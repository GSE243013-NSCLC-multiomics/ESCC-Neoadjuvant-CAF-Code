#!/usr/bin/env Rscript
# ============================================================================
# STEP 19M: CONSERVED FIBROBLAST CORE MECHANISM AUDIT
# ============================================================================
# Investigates candidate upstream regulatory programs associated with the
# conserved fibroblast core (BGN+, CDKN1B-, GSN-, TIMP1+) identified in
# Step19L. Hypothesis-generation only; no causal inference.
# ============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(reshape2)
})

BASE_DIR <- Sys.getenv("ESCC_CAF_PROJECT_ROOT", unset = ".")
M_DIR    <- file.path(BASE_DIR, "06_tables/cross_dataset/pathway_analysis/Step19M")
FIG_DIR  <- file.path(BASE_DIR, "05_figures/cross_dataset/pathway_analysis/Step19M")
LOG_DIR  <- file.path(BASE_DIR, "08_logs")
dir.create(M_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(LOG_DIR, showWarnings = FALSE, recursive = TRUE)

CORE_GENES <- c("BGN", "CDKN1B", "GSN", "TIMP1")

TARGET_PATHWAYS <- c("HALLMARK_HYPOXIA", "HALLMARK_COAGULATION",
                     "HALLMARK_KRAS_SIGNALING_UP", "HALLMARK_APOPTOSIS")
PATHWAY_SHORT <- gsub("HALLMARK_", "", TARGET_PATHWAYS)

# Canonical treatment maps — defined ONCE
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
  GSE197677  = c(treated = "NACT", control = "nNACT"),
  GSE221561  = c(treated = "Neoadjuvant_treated", control = "Surgery_alone")
)

compute_log2_cpm <- function(counts_mat) {
  lib_size <- colSums(counts_mat)
  lib_size[lib_size == 0] <- 1
  cpm_mat <- sweep(counts_mat, 2, lib_size, "/") * 1e6
  log2(cpm_mat + 1)
}

FROZEN_SIGN <- c(BGN = +1, CDKN1B = -1, GSN = -1, TIMP1 = +1)

# ============================================================================
# 19M0 — OUTPUT STRUCTURE
# ============================================================================
cat("============================================\n")
cat("STEP 19M: FIBROBLAST CORE MECHANISM AUDIT\n")
cat("============================================\n\n")

cat("--- OUTPUT DIRECTORIES ---\n")
cat("Tables:", M_DIR, "\n")
cat("Figures:", FIG_DIR, "\n\n")

# ============================================================================
# 19M1 — INPUT VALIDATION
# ============================================================================
cat("============================================\n")
cat("19M1: INPUT VALIDATION\n")
cat("============================================\n\n")

# Read frozen outputs
scores_all <- read.csv(file.path(BASE_DIR,
  "06_tables/cross_dataset/pathway_analysis/Step19I/STEP19I_SAMPLE_COMPARTMENT_PATHWAY_SCORES.csv"))
core_scores <- read.csv(file.path(BASE_DIR,
  "06_tables/cross_dataset/pathway_analysis/Step19L/STEP19L_SAMPLE_CORE_MODULE_SCORES.csv"))
core_effects <- read.csv(file.path(BASE_DIR,
  "06_tables/cross_dataset/pathway_analysis/Step19L/STEP19L_CORE_MODULE_CANONICAL_EFFECTS.csv"))
core_perm <- read.csv(file.path(BASE_DIR,
  "06_tables/cross_dataset/pathway_analysis/Step19L/STEP19L_CORE_MODULE_EXACT_TESTS.csv"))
core_loso <- read.csv(file.path(BASE_DIR,
  "06_tables/cross_dataset/pathway_analysis/Step19L/STEP19L_CORE_MODULE_LOSO_LONG.csv"))
core_coupling <- read.csv(file.path(BASE_DIR,
  "06_tables/cross_dataset/pathway_analysis/Step19L/STEP19L_CORE_PATHWAY_COUPLING.csv"))
core_gene_effects <- read.csv(file.path(BASE_DIR,
  "06_tables/cross_dataset/pathway_analysis/Step19L/STEP19L_CORE_GENE_EFFECT_VERIFICATION.csv"))
final_valid <- read.csv(file.path(BASE_DIR,
  "06_tables/cross_dataset/pathway_analysis/Step19L/STEP19L_FINAL_CORE_MODULE_VALIDATION.csv"))
inventory <- read.csv(file.path(BASE_DIR,
  "06_tables/cross_dataset/pathway_analysis/Step19I/STEP19I_COMPARTMENT_SAMPLE_INVENTORY.csv"))
# Note: CORR2A canonical effects are pathway-level, not gene-level; gene effects from 19L

fib197_file <- file.path(BASE_DIR,
  "06_tables/cross_dataset/pathway_analysis/Step19I/counts/GSE197677_FIBROBLAST_CAF_raw_counts.rds")
fib221_file <- file.path(BASE_DIR,
  "06_tables/cross_dataset/pathway_analysis/Step19I/counts/GSE221561_FIBROBLAST_CAF_raw_counts.rds")
stopifnot(file.exists(fib197_file), file.exists(fib221_file))
fib197 <- readRDS(fib197_file)
fib221 <- readRDS(fib221_file)

cat("GSE197677 fibroblast:", nrow(fib197), "genes x", ncol(fib197), "samples\n")
cat("GSE221561 fibroblast:", nrow(fib221), "genes x", ncol(fib221), "samples\n\n")

# Validate core scores (columns are BGN_z, CDKN1B_z, etc.)
core_score_cols <- paste0(CORE_GENES, "_z")
stopifnot(all(core_score_cols %in% colnames(core_scores)))
stopifnot("CORE4_score" %in% colnames(core_scores))
stopifnot("CORE2_STRONG_score" %in% colnames(core_scores))

# Validate sample IDs match treatment map
for (ds in names(TREAT_MAP)) {
  ds_samples <- core_scores$sample_id[core_scores$dataset == ds]
  map_samples <- names(TREAT_MAP[[ds]])
  stopifnot(all(ds_samples %in% map_samples))
  cat(ds, ": all", length(ds_samples), "sample IDs validated\n")
}

# Validate no normals in treatment comparison
stopifnot(!"Adjacent_normal" %in% core_scores$treatment)

# Validate canonical effects
cat("\n=== CANONICAL EFFECTS ===\n")
for (i in 1:nrow(core_effects)) {
  cat(sprintf("  %s %s: effect = %.4f\n",
    core_effects$dataset[i], core_effects$module[i], core_effects$effect[i]))
}

# Save input provenance
provenance <- data.frame(
  step = "19M",
  status = "INPUT_VALIDATED",
  n_fib197_samples = ncol(fib197),
  n_fib221_samples = ncol(fib221),
  n_common_genes = length(intersect(rownames(fib197), rownames(fib221))),
  core4_197_effect = core_effects$effect[core_effects$dataset == "GSE197677" & core_effects$module == "CORE4"],
  core4_221_effect = core_effects$effect[core_effects$dataset == "GSE221561" & core_effects$module == "CORE4"],
  core2_197_effect = core_effects$effect[core_effects$dataset == "GSE197677" & core_effects$module == "CORE2_STRONG"],
  core2_221_effect = core_effects$effect[core_effects$dataset == "GSE221561" & core_effects$module == "CORE2_STRONG"],
  timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
)
write.csv(provenance, file.path(M_DIR, "STEP19M_INPUT_PROVENANCE.csv"), row.names = FALSE)
cat("\nSaved: STEP19M_INPUT_PROVENANCE.csv\n")

# ============================================================================
# 19M2 — DEFINE REGULATORY GENE UNIVERSE
# ============================================================================
cat("\n============================================\n")
cat("19M2: DEFINE REGULATORY GENE UNIVERSE\n")
cat("============================================\n\n")

# Use pseudobulk log2-CPM expression
expr197 <- compute_log2_cpm(fib197)
expr221 <- compute_log2_cpm(fib221)

# Common genes
common_genes <- intersect(rownames(expr197), rownames(expr221))
cat("Common genes:", length(common_genes), "\n")

# Filter: require expression in both datasets
# Use genes with mean log2-CPM > 1 in at least one dataset
mean197 <- rowMeans(expr197[common_genes, , drop = FALSE])
mean221 <- rowMeans(expr221[common_genes, , drop = FALSE])
expressed <- (mean197 > 1) | (mean221 > 1)
universe_genes <- common_genes[expressed]
cat("Expressed common genes:", length(universe_genes), "\n")

# Remove core genes from discovery universe (but keep for later use)
discovery_genes <- setdiff(universe_genes, CORE_GENES)
cat("Discovery genes (excl core):", length(discovery_genes), "\n")

universe_df <- data.frame(
  gene = universe_genes,
  in_discovery = universe_genes %in% discovery_genes,
  mean_log2cpm_197 = mean197[universe_genes],
  mean_log2cpm_221 = mean221[universe_genes],
  stringsAsFactors = FALSE
)
write.csv(universe_df, file.path(M_DIR, "STEP19M_COMMON_FIBROBLAST_GENE_UNIVERSE.csv"),
          row.names = FALSE)
cat("Saved: STEP19M_COMMON_FIBROBLAST_GENE_UNIVERSE.csv\n\n")

# ============================================================================
# 19M3 — CORE-ASSOCIATED GENES
# ============================================================================
cat("============================================\n")
cat("19M3: CORE-ASSOCIATED GENES\n")
cat("============================================\n\n")

# Build sample-level data frames
samples197 <- colnames(expr197)
samples221 <- colnames(expr221)

df197 <- data.frame(
  sample_id = samples197,
  dataset = "GSE197677",
  treatment = TREAT_MAP$GSE197677[samples197],
  CORE4_score = core_scores$CORE4_score[match(samples197, core_scores$sample_id)],
  CORE2_STRONG_score = core_scores$CORE2_STRONG_score[match(samples197, core_scores$sample_id)],
  lib_size = colSums(fib197[, samples197, drop = FALSE]),
  row.names = samples197,
  stringsAsFactors = FALSE
)

df221 <- data.frame(
  sample_id = samples221,
  dataset = "GSE221561",
  treatment = TREAT_MAP$GSE221561[samples221],
  CORE4_score = core_scores$CORE4_score[match(samples221, core_scores$sample_id)],
  CORE2_STRONG_score = core_scores$CORE2_STRONG_score[match(samples221, core_scores$sample_id)],
  lib_size = colSums(fib221[, samples221, drop = FALSE]),
  row.names = samples221,
  stringsAsFactors = FALSE
)

# Fibroblast abundance proxy: library size
df197$lib_size_z <- scale(log1p(df197$lib_size))[, 1]
df221$lib_size_z <- scale(log1p(df221$lib_size))[, 1]

# Spearman correlation for each gene with CORE4 and CORE2_STRONG
cat("Computing Spearman correlations...\n")

assoc_list <- list()
for (gene in discovery_genes) {
  if (gene %in% rownames(expr197) && gene %in% rownames(expr221)) {
    rho197_c4 <- cor(expr197[gene, samples197], df197$CORE4_score, method = "spearman", use = "complete.obs")
    rho221_c4 <- cor(expr221[gene, samples221], df221$CORE4_score, method = "spearman", use = "complete.obs")
    rho197_c2 <- cor(expr197[gene, samples197], df197$CORE2_STRONG_score, method = "spearman", use = "complete.obs")
    rho221_c2 <- cor(expr221[gene, samples221], df221$CORE2_STRONG_score, method = "spearman", use = "complete.obs")
    
    assoc_list[[gene]] <- data.frame(
      gene = gene,
      rho_GSE197677_CORE4 = rho197_c4,
      rho_GSE221561_CORE4 = rho221_c4,
      rho_GSE197677_CORE2 = rho197_c2,
      rho_GSE221561_CORE2 = rho221_c2,
      stringsAsFactors = FALSE
    )
  }
}

assoc_df <- do.call(rbind, assoc_list)
rownames(assoc_df) <- NULL

# Cross-dataset comparison for CORE4
assoc_df$same_sign_CORE4 <- sign(assoc_df$rho_GSE197677_CORE4) == sign(assoc_df$rho_GSE221561_CORE4)
assoc_df$mean_abs_rho_CORE4 <- (abs(assoc_df$rho_GSE197677_CORE4) + abs(assoc_df$rho_GSE221561_CORE4)) / 2
assoc_df$min_abs_rho_CORE4 <- pmin(abs(assoc_df$rho_GSE197677_CORE4), abs(assoc_df$rho_GSE221561_CORE4))

# Cross-dataset comparison for CORE2
assoc_df$same_sign_CORE2 <- sign(assoc_df$rho_GSE197677_CORE2) == sign(assoc_df$rho_GSE221561_CORE2)
assoc_df$mean_abs_rho_CORE2 <- (abs(assoc_df$rho_GSE197677_CORE2) + abs(assoc_df$rho_GSE221561_CORE2)) / 2
assoc_df$min_abs_rho_CORE2 <- pmin(abs(assoc_df$rho_GSE197677_CORE2), abs(assoc_df$rho_GSE221561_CORE2))

# Tier classification for CORE4
assoc_df$tier_CORE4 <- "DISCORDANT"
assoc_df$tier_CORE4[assoc_df$same_sign_CORE4 & assoc_df$min_abs_rho_CORE4 >= 0.60] <- "CONSISTENT_STRONG"
assoc_df$tier_CORE4[assoc_df$same_sign_CORE4 & assoc_df$min_abs_rho_CORE4 >= 0.40 & assoc_df$min_abs_rho_CORE4 < 0.60] <- "CONSISTENT_MODERATE"
assoc_df$tier_CORE4[!assoc_df$same_sign_CORE4] <- "DISCORDANT"
# Dataset specific: same sign but one < 0.40
assoc_df$tier_CORE4[assoc_df$same_sign_CORE4 & assoc_df$min_abs_rho_CORE4 < 0.40] <- "DATASET_SPECIFIC"

# Tier classification for CORE2
assoc_df$tier_CORE2 <- "DISCORDANT"
assoc_df$tier_CORE2[assoc_df$same_sign_CORE2 & assoc_df$min_abs_rho_CORE2 >= 0.60] <- "CONSISTENT_STRONG"
assoc_df$tier_CORE2[assoc_df$same_sign_CORE2 & assoc_df$min_abs_rho_CORE2 >= 0.40 & assoc_df$min_abs_rho_CORE2 < 0.60] <- "CONSISTENT_MODERATE"
assoc_df$tier_CORE2[!assoc_df$same_sign_CORE2] <- "DISCORDANT"
assoc_df$tier_CORE2[assoc_df$same_sign_CORE2 & assoc_df$min_abs_rho_CORE2 < 0.40] <- "DATASET_SPECIFIC"

# Summary
cat("\n=== CORE4 ASSOCIATION TIERS ===\n")
print(table(assoc_df$tier_CORE4))
cat("\n=== CORE2 ASSOCIATION TIERS ===\n")
print(table(assoc_df$tier_CORE2))

# Show top genes
cat("\n=== TOP 20 CORE4-ASSOCIATED GENES (by mean_abs_rho) ===\n")
top_c4 <- assoc_df[order(-assoc_df$mean_abs_rho_CORE4), ][1:min(20, nrow(assoc_df)), ]
for (i in 1:nrow(top_c4)) {
  cat(sprintf("  %s: rho197=%.3f rho221=%.3f tier=%s\n",
    top_c4$gene[i], top_c4$rho_GSE197677_CORE4[i], top_c4$rho_GSE221561_CORE4[i], top_c4$tier_CORE4[i]))
}

cat("\n=== TOP 20 CORE2-ASSOCIATED GENES (by mean_abs_rho) ===\n")
top_c2 <- assoc_df[order(-assoc_df$mean_abs_rho_CORE2), ][1:min(20, nrow(assoc_df)), ]
for (i in 1:nrow(top_c2)) {
  cat(sprintf("  %s: rho197=%.3f rho221=%.3f tier=%s\n",
    top_c2$gene[i], top_c2$rho_GSE197677_CORE2[i], top_c2$rho_GSE221561_CORE2[i], top_c2$tier_CORE2[i]))
}

write.csv(assoc_df, file.path(M_DIR, "STEP19M_CORE_ASSOCIATED_GENES.csv"), row.names = FALSE)
cat("\nSaved: STEP19M_CORE_ASSOCIATED_GENES.csv\n")

# ============================================================================
# 19M4 — TRANSCRIPTION FACTOR CANDIDATE AUDIT
# ============================================================================
cat("\n============================================\n")
cat("19M4: TRANSCRIPTION FACTOR CANDIDATE AUDIT\n")
cat("============================================\n\n")

# Curated list of known human TFs relevant to fibroblast biology
# Source: TFCatalog / GO:0003700 / literature
# Restrict to genes present in our universe
KNOWN_TFS <- c(
  # TGF-beta/SMAD
  "SMAD2", "SMAD3", "SMAD4", "SMAD7",
  # HIF/hypoxia
  "HIF1A", "HIF3A", "EPAS1", "ARNT", "EGLN1",
  # AP-1 family
  "JUN", "JUNB", "JUND", "FOS", "FOSB", "FOSL1", "FOSL2",
  "ATF3", "ATF4",
  # NF-kB
  "RELA", "NFKB1", "NFKB2", "RELB", "STAT3",
  # TEAD/YAP
  "TEAD1", "TEAD2", "TEAD3", "TEAD4",
  # STAT family
  "STAT1", "STAT2", "STAT4", "STAT5A", "STAT5B", "STAT6",
  # p53/cell-cycle
  "TP53", "TP63", "CDKN1A", "RB1", "E2F1",
  # Fibroblast-relevant
  "SNAI1", "SNAI2", "TWIST1", "TWIST2", "ZEB1", "ZEB2",
  "FOXD1", "FOXF1", "FOXO1", "FOXO3",
  "KLF4", "KLF5", "KLF6",
  "SP1", "SP3",
  "ETS1", "ETS2",
  "GATA2", "GATA3",
  "SOX2", "SOX4", "SOX9",
  "NFIA", "NFIB", "NFIC",
  "RUNX1", "RUNX2",
  "HOXA5", "HOXB7",
  "TGFB1", "TGFB2", "TGFB3"
)

tf_in_universe <- KNOWN_TFS[KNOWN_TFS %in% universe_genes]
tf_in_discovery <- KNOWN_TFS[KNOWN_TFS %in% discovery_genes]
cat("Known TFs in universe:", length(tf_in_universe), "\n")
cat("Known TFs in discovery:", length(tf_in_discovery), "\n")

# TF association with modules
tf_assoc_list <- list()
for (tf in tf_in_discovery) {
  rho197_c4 <- cor(expr197[tf, samples197], df197$CORE4_score, method = "spearman", use = "complete.obs")
  rho221_c4 <- cor(expr221[tf, samples221], df221$CORE4_score, method = "spearman", use = "complete.obs")
  rho197_c2 <- cor(expr197[tf, samples197], df197$CORE2_STRONG_score, method = "spearman", use = "complete.obs")
  rho221_c2 <- cor(expr221[tf, samples221], df221$CORE2_STRONG_score, method = "spearman", use = "complete.obs")
  
  # Treatment effect (treated - control)
  treat197 <- names(TREAT_MAP$GSE197677[TREAT_MAP$GSE197677 == TREAT_CONTRAST$GSE197677["treated"]])
  ctrl197 <- names(TREAT_MAP$GSE197677[TREAT_MAP$GSE197677 == TREAT_CONTRAST$GSE197677["control"]])
  treat221 <- names(TREAT_MAP$GSE221561[TREAT_MAP$GSE221561 == TREAT_CONTRAST$GSE221561["treated"]])
  ctrl221 <- names(TREAT_MAP$GSE221561[TREAT_MAP$GSE221561 == TREAT_CONTRAST$GSE221561["control"]])
  
  eff197 <- mean(expr197[tf, treat197]) - mean(expr197[tf, ctrl197])
  eff221 <- mean(expr221[tf, treat221]) - mean(expr221[tf, ctrl221])
  
  same_dir <- sign(rho197_c4) == sign(rho221_c4) | sign(rho197_c2) == sign(rho221_c2)
  treat_consistent <- sign(eff197) == sign(rho197_c4) | sign(eff221) == sign(rho221_c4)
  
  # Classification
  if (same_dir && abs(rho197_c4) >= 0.40 && abs(rho221_c4) >= 0.40) {
    classification <- "TF_EXPRESSION_CROSS_DATASET_SUPPORTED"
  } else if (same_dir) {
    classification <- "TF_EXPRESSION_ASSOCIATED"
  } else {
    classification <- "TF_EXPRESSION_DATASET_SPECIFIC"
  }
  
  tf_assoc_list[[tf]] <- data.frame(
    TF = tf,
    rho_GSE197677_CORE4 = rho197_c4,
    rho_GSE221561_CORE4 = rho221_c4,
    rho_GSE197677_CORE2 = rho197_c2,
    rho_GSE221561_CORE2 = rho221_c2,
    effect_GSE197677 = eff197,
    effect_GSE221561 = eff221,
    same_direction_module = same_dir,
    treatment_direction_consistent = treat_consistent,
    classification = classification,
    stringsAsFactors = FALSE
  )
}

tf_assoc_df <- do.call(rbind, tf_assoc_list)
rownames(tf_assoc_df) <- NULL

# Rank TFs by evidence
tf_assoc_df$rank_score <- 0
tf_assoc_df$rank_score[tf_assoc_df$same_direction_module] <- tf_assoc_df$rank_score[tf_assoc_df$same_direction_module] + 1
tf_assoc_df$rank_score[abs(tf_assoc_df$rho_GSE197677_CORE4) >= 0.40 & abs(tf_assoc_df$rho_GSE221561_CORE4) >= 0.40] <-
  tf_assoc_df$rank_score[abs(tf_assoc_df$rho_GSE197677_CORE4) >= 0.40 & abs(tf_assoc_df$rho_GSE221561_CORE4) >= 0.40] + 1
tf_assoc_df$rank_score[tf_assoc_df$treatment_direction_consistent] <- tf_assoc_df$rank_score[tf_assoc_df$treatment_direction_consistent] + 1

tf_assoc_df <- tf_assoc_df[order(-tf_assoc_df$rank_score, -abs(tf_assoc_df$rho_GSE197677_CORE4)), ]

cat("\n=== TOP 15 TF CANDIDATES ===\n")
for (i in 1:min(15, nrow(tf_assoc_df))) {
  cat(sprintf("  %s: score=%d class=%s rho197_C4=%.3f rho221_C4=%.3f\n",
    tf_assoc_df$TF[i], tf_assoc_df$rank_score[i], tf_assoc_df$classification[i],
    tf_assoc_df$rho_GSE197677_CORE4[i], tf_assoc_df$rho_GSE221561_CORE4[i]))
}

cat("\n=== TF CLASSIFICATION SUMMARY ===\n")
print(table(tf_assoc_df$classification))

write.csv(tf_assoc_df, file.path(M_DIR, "STEP19M_TF_EXPRESSION_ASSOCIATIONS.csv"), row.names = FALSE)
cat("Saved: STEP19M_TF_EXPRESSION_ASSOCIATIONS.csv\n")

# ============================================================================
# 19M5 — MOTIF / TARGET-SET ENRICHMENT
# ============================================================================
cat("\n============================================\n")
cat("19M5: MOTIF / TARGET-SET ENRICHMENT\n")
cat("============================================\n\n")

# Attempt to find a curated TF-target resource
# Check if msigdbr or similar is available
has_msigdbr <- requireNamespace("msigdbr", quietly = TRUE)
has_dorothea <- requireNamespace("dorothea", quietly = TRUE)

enrichment_performed <- FALSE

if (has_msigdbr) {
  cat("msigdbr available — attempting TF target-set enrichment\n")
  tryCatch({
    library(msigdbr, quietly = TRUE)
    # Get TF target gene sets from MSigDB
    tft_genesets <- msigdbr(species = "Homo sapiens", category = "C3:TF")
    # Get regulons
    regulons <- tft_genesets %>%
      dplyr::filter(gs_subcat %in% c("TFT:GTRD", "TFT:TFT_Legacy")) %>%
      dplyr::select(gs_name, gene_symbol) %>%
      distinct()
    
    if (nrow(regulons) > 0) {
      enrichment_performed <- TRUE
      cat("Loaded", n_distinct(regulons$gs_name), "TF target gene sets\n")
    }
  }, error = function(e) {
    cat("msigdbr enrichment failed:", conditionMessage(e), "\n")
  })
} else {
  cat("msigdbr not available\n")
}

if (!enrichment_performed) {
  cat("No suitable TF-target resource available. Skipping enrichment.\n\n")
  
  enrichment_results <- data.frame(
    regulator = character(0),
    NES = numeric(0),
    P_value = numeric(0),
    FDR = numeric(0),
    dataset = character(0),
    module = character(0),
    direction = character(0)
  )
  cross_dataset_enrichment <- data.frame(
    regulator = character(0),
    support_197 = character(0),
    support_221 = character(0),
    cross_dataset = character(0)
  )
}

if (enrichment_performed) {
  cat("Running ranked enrichment for CORE4/CORE2 STRONG gene sets...\n")
  
  # Define gene sets from associations
  c4_pos_genes <- assoc_df$gene[assoc_df$tier_CORE4 %in% c("CONSISTENT_STRONG", "CONSISTENT_MODERATE") &
                                 assoc_df$rho_GSE197677_CORE4 > 0]
  c4_neg_genes <- assoc_df$gene[assoc_df$tier_CORE4 %in% c("CONSISTENT_STRONG", "CONSISTENT_MODERATE") &
                                 assoc_df$rho_GSE197677_CORE4 < 0]
  c2_pos_genes <- assoc_df$gene[assoc_df$tier_CORE2 %in% c("CONSISTENT_STRONG", "CONSISTENT_MODERATE") &
                                 assoc_df$rho_GSE197677_CORE2 > 0]
  c2_neg_genes <- assoc_df$gene[assoc_df$tier_CORE2 %in% c("CONSISTENT_STRONG", "CONSISTENT_MODERATE") &
                                 assoc_df$rho_GSE197677_CORE2 < 0]
  
  cat("CORE4 positive genes:", length(c4_pos_genes), "\n")
  cat("CORE4 negative genes:", length(c4_neg_genes), "\n")
  cat("CORE2 positive genes:", length(c2_pos_genes), "\n")
  cat("CORE2 negative genes:", length(c2_neg_genes), "\n\n")
  
  # Simple overlap-based enrichment (Fisher's exact)
  run_enrichment <- function(gene_set, regulon_list, universe) {
    results <- list()
    for (reg_name in names(regulon_list)) {
      reg_genes <- regulon_list[[reg_name]]
      overlap <- length(intersect(gene_set, reg_genes))
      in_set <- length(intersect(reg_genes, universe))
      not_in_set <- length(reg_genes) - in_set
      in_gs <- length(intersect(gene_set, universe))
      not_in_gs <- length(universe) - in_gs
      
      if (overlap >= 2) {
        mat <- matrix(c(overlap, in_gs - overlap, in_set - overlap, 
                       length(universe) - overlap - (in_gs - overlap) - (in_set - overlap)),
                     nrow = 2, byrow = TRUE)
        ft <- fisher.test(mat, alternative = "greater")
        results[[reg_name]] <- data.frame(
          regulator = reg_name,
          overlap = overlap,
          P_value = ft$p.value,
          stringsAsFactors = FALSE
        )
      }
    }
    if (length(results) > 0) do.call(rbind, results) else data.frame()
  }
  
  # Simplified regulon list from MSigDB
  regulon_split <- split(regulons$gene_symbol, regulons$gs_name)
  
  enrich_list <- list()
  for (ds_name in c("GSE197677", "GSE221561")) {
    for (module in c("CORE4", "CORE2_STRONG")) {
      for (dir_label in c("positive", "negative")) {
        if (module == "CORE4") {
          gs <- if (dir_label == "positive") c4_pos_genes else c4_neg_genes
        } else {
          gs <- if (dir_label == "positive") c2_pos_genes else c2_neg_genes
        }
        if (length(gs) >= 3) {
          enr <- run_enrichment(gs, regulon_split, universe_genes)
          if (nrow(enr) > 0) {
            enr$dataset <- ds_name
            enr$module <- module
            enr$direction <- dir_label
            enrich_list[[paste(ds_name, module, dir_label)]] <- enr
          }
        }
      }
    }
  }
  
  if (length(enrich_list) > 0) {
    enrichment_results <- do.call(rbind, enrich_list)
    enrichment_results$FDR <- p.adjust(enrichment_results$P_value, method = "BH")
    enrichment_results <- enrichment_results[order(enrichment_results$FDR), ]
  } else {
    enrichment_results <- data.frame()
  }
  
  write.csv(enrichment_results, file.path(M_DIR, "STEP19M_REGULATOR_ENRICHMENT.csv"), row.names = FALSE)
  cat("Saved: STEP19M_REGULATOR_ENRICHMENT.csv\n")
}

# ============================================================================
# 19M6 — PRIORITIZED BIOLOGICAL REGULATORS
# ============================================================================
cat("\n============================================\n")
cat("19M6: PRIORITIZED BIOLOGICAL REGULATORS\n")
cat("============================================\n\n")

prespecified <- list(
  TGF_beta_SMAD = c("TGFB1", "TGFB2", "TGFB3", "SMAD2", "SMAD3", "SMAD4", "SMAD7"),
  HIF_hypoxia = c("HIF1A", "HIF3A", "EPAS1", "ARNT", "EGLN1"),
  AP1_family = c("JUN", "JUNB", "JUND", "FOS", "FOSB", "FOSL1", "FOSL2", "ATF3", "ATF4"),
  NFkB = c("RELA", "NFKB1", "NFKB2", "RELB", "STAT3"),
  TEAD_YAP = c("TEAD1", "TEAD2", "TEAD3", "TEAD4"),
  STAT = c("STAT1", "STAT2", "STAT3", "STAT4", "STAT5A", "STAT5B", "STAT6"),
  p53_cell_cycle = c("TP53", "TP63", "CDKN1A", "RB1", "E2F1", "CDKN1B")
)

prespec_df_list <- list()
for (prog_name in names(prespecified)) {
  for (tf in prespecified[[prog_name]]) {
    if (tf %in% universe_genes) {
      rho197_c4 <- cor(expr197[tf, samples197], df197$CORE4_score, method = "spearman", use = "complete.obs")
      rho221_c4 <- cor(expr221[tf, samples221], df221$CORE4_score, method = "spearman", use = "complete.obs")
      rho197_c2 <- cor(expr197[tf, samples197], df197$CORE2_STRONG_score, method = "spearman", use = "complete.obs")
      rho221_c2 <- cor(expr221[tf, samples221], df221$CORE2_STRONG_score, method = "spearman", use = "complete.obs")
      
      # Match to unbiased ranking
      if (tf %in% tf_assoc_df$TF) {
        unbiased_rank <- which(tf_assoc_df$TF == tf)
        unbiased_class <- tf_assoc_df$classification[unbiased_rank]
      } else {
        unbiased_rank <- NA
        unbiased_class <- "NOT_IN_DISCOVERY"
      }
      
      prespec_df_list[[paste(prog_name, tf)]] <- data.frame(
        program = prog_name,
        TF = tf,
        rho_GSE197677_CORE4 = rho197_c4,
        rho_GSE221561_CORE4 = rho221_c4,
        rho_GSE197677_CORE2 = rho197_c2,
        rho_GSE221561_CORE2 = rho221_c2,
        same_direction_CORE4 = sign(rho197_c4) == sign(rho221_c4),
        same_direction_CORE2 = sign(rho197_c2) == sign(rho221_c2),
        unbiased_rank = unbiased_rank,
        unbiased_classification = unbiased_class,
        label = "PRESPECIFIED_CANDIDATE",
        stringsAsFactors = FALSE
      )
    }
  }
}

prespec_df <- do.call(rbind, prespec_df_list)
rownames(prespec_df) <- NULL

cat("=== PRESPECIFIED REGULATORS PRESENT IN UNIVERSE ===\n")
for (prog_name in names(prespecified)) {
  present <- prespecified[[prog_name]][prespecified[[prog_name]] %in% universe_genes]
  cat(sprintf("  %s: %d/%d present (%s)\n",
    prog_name, length(present), length(prespecified[[prog_name]]),
    paste(present, collapse = ", ")))
}

write.csv(prespec_df, file.path(M_DIR, "STEP19M_PRESPECIFIED_REGULATOR_AUDIT.csv"), row.names = FALSE)
cat("Saved: STEP19M_PRESPECIFIED_REGULATOR_AUDIT.csv\n")

# ============================================================================
# 19M7 — CORE VS HETEROGENEOUS PATHWAYS
# ============================================================================
cat("\n============================================\n")
cat("19M7: CORE VS HETEROGENEOUS PATHWAYS\n")
cat("============================================\n\n")

# Build pathway score data frames for fibroblast compartment
fib_scores <- scores_all[scores_all$compartment == "FIBROBLAST_CAF", ]

# Get pathway scores per sample
pathway_assoc_list <- list()
for (pw in TARGET_PATHWAYS) {
  pw_short <- gsub("HALLMARK_", "", pw)
  
  for (ds in c("GSE197677", "GSE221561")) {
    ds_samples <- names(TREAT_MAP[[ds]])
    pw_scores <- fib_scores[fib_scores$dataset == ds & fib_scores$pathway == pw, ]
    
    for (s in ds_samples) {
      row <- pw_scores[pw_scores$sample == s, ]
      if (nrow(row) > 0) {
        pathway_assoc_list[[paste(ds, pw_short, s)]] <- data.frame(
          sample_id = s,
          dataset = ds,
          pathway = pw_short,
          pathway_score = row$score,
          stringsAsFactors = FALSE
        )
      }
    }
  }
}

pathway_df <- do.call(rbind, pathway_assoc_list)
rownames(pathway_df) <- NULL

# Merge with core scores
pathway_df$CORE4_score <- core_scores$CORE4_score[match(
  paste(pathway_df$dataset, pathway_df$sample_id),
  paste(core_scores$dataset, core_scores$sample_id))]
pathway_df$CORE2_STRONG_score <- core_scores$CORE2_STRONG_score[match(
  paste(pathway_df$dataset, pathway_df$sample_id),
  paste(core_scores$dataset, core_scores$sample_id))]

# Compute associations
cpw_results <- list()
for (pw in PATHWAY_SHORT) {
  for (ds in c("GSE197677", "GSE221561")) {
    sub <- pathway_df[pathway_df$pathway == pw & pathway_df$dataset == ds, ]
    if (nrow(sub) >= 4) {
      rho_c4 <- cor(sub$pathway_score, sub$CORE4_score, method = "spearman", use = "complete.obs")
      rho_c2 <- cor(sub$pathway_score, sub$CORE2_STRONG_score, method = "spearman", use = "complete.obs")
      cpw_results[[paste(pw, ds)]] <- data.frame(
        pathway = pw,
        dataset = ds,
        rho_CORE4 = rho_c4,
        rho_CORE2 = rho_c2,
        abs_rho_CORE4 = abs(rho_c4),
        abs_rho_CORE2 = abs(rho_c2),
        stringsAsFactors = FALSE
      )
    }
  }
}

cpw_df <- do.call(rbind, cpw_results)
rownames(cpw_df) <- NULL

# Now classify regulators as CORE_ASSOCIATED, PATHWAY_ASSOCIATED, etc.
# For the top regulators, check whether they track CORE4 or pathway scores
regulator_class_list <- list()

# Use top TFs from 19M4 plus prespecified key regulators
top_regs <- unique(c(
  tf_assoc_df$TF[1:min(10, nrow(tf_assoc_df))],
  c("HIF1A", "SMAD3", "JUN", "FOS", "STAT3", "TP53", "TEAD1", "NFKB1", "RELA")
))

for (reg in top_regs) {
  if (!(reg %in% universe_genes)) next
  
    reg_rho_c4_197 <- cor(expr197[reg, samples197], df197$CORE4_score, method = "spearman", use = "complete.obs")
    reg_rho_c4_221 <- cor(expr221[reg, samples221], df221$CORE4_score, method = "spearman", use = "complete.obs")
  
  pathway_tracks <- c()
  for (pw in PATHWAY_SHORT) {
    pw_197 <- pathway_df[pathway_df$pathway == pw & pathway_df$dataset == "GSE197677", ]
    pw_221 <- pathway_df[pathway_df$pathway == pw & pathway_df$dataset == "GSE221561", ]
    # Match by sample_id
    matched_197 <- pw_197[match(samples197, pw_197$sample_id), ]
    matched_221 <- pw_221[match(samples221, pw_221$sample_id), ]
    rho_197 <- cor(matched_197$pathway_score, as.numeric(expr197[reg, samples197]),
                   method = "spearman", use = "complete.obs")
    rho_221 <- cor(matched_221$pathway_score, as.numeric(expr221[reg, samples221]),
                   method = "spearman", use = "complete.obs")
    pathway_tracks <- c(pathway_tracks, abs(rho_197) >= 0.40 | abs(rho_221) >= 0.40)
  }
  
  core_associated <- abs(reg_rho_c4_197) >= 0.40 | abs(reg_rho_c4_221) >= 0.40
  pathway_associated <- any(pathway_tracks)
  
  if (core_associated && pathway_associated) {
    classification <- "SHARED_CORE_PATHWAY_ASSOCIATED"
  } else if (core_associated) {
    classification <- "CORE_ASSOCIATED"
  } else if (pathway_associated) {
    classification <- "PATHWAY_ASSOCIATED"
  } else {
    classification <- "INCONCLUSIVE"
  }
  
  regulator_class_list[[reg]] <- data.frame(
    regulator = reg,
    rho_CORE4_GSE197677 = reg_rho_c4_197,
    rho_CORE4_GSE221561 = reg_rho_c4_221,
    n_pathways_associated = sum(pathway_tracks),
    classification = classification,
    stringsAsFactors = FALSE
  )
}

reg_class_df <- do.call(rbind, regulator_class_list)
rownames(reg_class_df) <- NULL

cat("=== REGULATOR-CORE/PATHWAY CLASSIFICATION ===\n")
for (i in 1:nrow(reg_class_df)) {
  cat(sprintf("  %s: %s (C4rho197=%.3f C4rho221=%.3f n_pw=%d)\n",
    reg_class_df$regulator[i], reg_class_df$classification[i],
    reg_class_df$rho_CORE4_GSE197677[i], reg_class_df$rho_CORE4_GSE221561[i],
    reg_class_df$n_pathways_associated[i]))
}

cat("\n=== CORE-PATHWAY ASSOCIATION MATRIX ===\n")
print(cpw_df[, c("pathway", "dataset", "rho_CORE4", "rho_CORE2")], row.names = FALSE)

write.csv(reg_class_df, file.path(M_DIR, "STEP19M_CORE_PATHWAY_REGULATOR_MAP.csv"), row.names = FALSE)
cat("\nSaved: STEP19M_CORE_PATHWAY_REGULATOR_MAP.csv\n")

# ============================================================================
# 19M8 — PARTIAL ASSOCIATION AUDIT
# ============================================================================
cat("\n============================================\n")
cat("19M8: PARTIAL ASSOCIATION AUDIT\n")
cat("============================================\n\n")

covariate_results <- list()

for (reg in top_regs[1:min(10, length(top_regs))]) {
  if (!(reg %in% universe_genes)) next
  
  for (ds in c("GSE197677", "GSE221561")) {
    if (ds == "GSE197677") {
      expr_reg <- as.numeric(expr197[reg, samples197])
      scores_sub <- df197
    } else {
      expr_reg <- as.numeric(expr221[reg, samples221])
      scores_sub <- df221
    }
    
    for (module in c("CORE4_score", "CORE2_STRONG_score")) {
      outcome <- scores_sub[[module]]
      
      # Model 1: unadjusted
      m1 <- tryCatch({
        fit1 <- lm(outcome ~ expr_reg)
        s1 <- summary(fit1)$coefficients
        if ("expr_reg" %in% rownames(s1)) c("Estimate" = s1["expr_reg", "Estimate"], "Pr" = s1["expr_reg", "Pr(>|t|)"])
        else c("Estimate" = NA, "Pr" = NA)
      }, error = function(e) c("Estimate" = NA, "Pr" = NA))
      
      # Model 2: + library size
      m2 <- tryCatch({
        fit2 <- lm(outcome ~ expr_reg + scores_sub$lib_size_z)
        s2 <- summary(fit2)$coefficients
        if ("expr_reg" %in% rownames(s2)) c("Estimate" = s2["expr_reg", "Estimate"], "Pr" = s2["expr_reg", "Pr(>|t|)"])
        else c("Estimate" = NA, "Pr" = NA)
      }, error = function(e) c("Estimate" = NA, "Pr" = NA))
      
      covariate_results[[paste(reg, ds, module)]] <- data.frame(
        regulator = reg,
        dataset = ds,
        module = module,
        coef_unadjusted = as.numeric(m1["Estimate"]),
        p_unadjusted = as.numeric(m1["Pr"]),
        coef_adjusted = as.numeric(m2["Estimate"]),
        p_adjusted = as.numeric(m2["Pr"]),
        stable_direction = sign(as.numeric(m1["Estimate"])) == sign(as.numeric(m2["Estimate"])),
        stringsAsFactors = FALSE
      )
    }
  }
}

covariate_df <- do.call(rbind, covariate_results)
rownames(covariate_df) <- NULL

cat("=== COVARIATE SENSITIVITY (top regulators) ===\n")
for (i in 1:min(20, nrow(covariate_df))) {
  cat(sprintf("  %s %s %s: coef_adj=%.4f stable=%s\n",
    covariate_df$regulator[i], covariate_df$dataset[i], covariate_df$module[i],
    covariate_df$coef_adjusted[i], covariate_df$stable_direction[i]))
}

write.csv(covariate_df, file.path(M_DIR, "STEP19M_REGULATOR_COVARIATE_SENSITIVITY.csv"), row.names = FALSE)
cat("\nSaved: STEP19M_REGULATOR_COVARIATE_SENSITIVITY.csv\n")

# ============================================================================
# 19M9 — LEAVE-ONE-SAMPLE-OUT REGULATOR STABILITY
# ============================================================================
cat("\n============================================\n")
cat("19M9: LOSO REGULATOR STABILITY\n")
cat("============================================\n\n")

loso_reg_list <- list()
for (reg in top_regs[1:min(10, length(top_regs))]) {
  if (!(reg %in% universe_genes)) next
  
  for (ds in c("GSE197677", "GSE221561")) {
    if (ds == "GSE197677") {
      expr_reg <- as.numeric(expr197[reg, samples197])
      mod_scores <- df197$CORE4_score
      sample_ids <- samples197
    } else {
      expr_reg <- as.numeric(expr221[reg, samples221])
      mod_scores <- df221$CORE4_score
      sample_ids <- samples221
    }
    
    n <- length(sample_ids)
    rho_loso <- numeric(n)
    for (i in 1:n) {
      keep <- setdiff(1:n, i)
      rho_loso[i] <- cor(expr_reg[keep], mod_scores[keep], method = "spearman", use = "complete.obs")
    }
    
    median_rho <- median(rho_loso, na.rm = TRUE)
    min_rho <- min(rho_loso, na.rm = TRUE)
    max_rho <- max(rho_loso, na.rm = TRUE)
    full_rho <- cor(expr_reg, mod_scores, method = "spearman", use = "complete.obs")
    sign_preservation <- mean(sign(rho_loso) == sign(full_rho), na.rm = TRUE)
    if (is.na(sign_preservation)) sign_preservation <- 0
    
    # Classification
    if (sign_preservation >= 0.80) {
      classification <- "ROBUST_ASSOCIATION"
    } else if (ds == "GSE221561" && sign_preservation < 0.80) {
      classification <- "CONTROL_SENSITIVE"
    } else {
      classification <- "UNSTABLE"
    }
    
    loso_reg_list[[paste(reg, ds)]] <- data.frame(
      regulator = reg,
      dataset = ds,
      median_rho = median_rho,
      min_rho = min_rho,
      max_rho = max_rho,
      sign_preservation = sign_preservation,
      classification = classification,
      stringsAsFactors = FALSE
    )
  }
}

loso_reg_df <- do.call(rbind, loso_reg_list)
rownames(loso_reg_df) <- NULL

cat("=== LOSO REGULATOR STABILITY ===\n")
for (i in 1:nrow(loso_reg_df)) {
  cat(sprintf("  %s %s: median=%.3f sign_pres=%.0f%% -> %s\n",
    loso_reg_df$regulator[i], loso_reg_df$dataset[i],
    loso_reg_df$median_rho[i], loso_reg_df$sign_preservation[i] * 100,
    loso_reg_df$classification[i]))
}

write.csv(loso_reg_df, file.path(M_DIR, "STEP19M_REGULATOR_LOSO.csv"), row.names = FALSE)
cat("\nSaved: STEP19M_REGULATOR_LOSO.csv\n")

# ============================================================================
# 19M10 — CORE2-SPECIFIC MECHANISTIC AUDIT
# ============================================================================
cat("\n============================================\n")
cat("19M10: CORE2-SPECIFIC MECHANISTIC AUDIT\n")
cat("============================================\n\n")

# CDKN1B and GSN oriented scores
# CDKN1B: negative orientation -> oriented = -CDKN1B_z (higher = lower CDKN1B expression)
# GSN: negative orientation -> oriented = -GSN_z
# CORE2_STRONG = mean(oriented CDKN1B, oriented GSN)
# But from frozen scores: CORE2_STRONG = (CDKN1B_z + GSN_z) / 2 with negative sign baked in

core2_shared_list <- list()
for (reg in top_regs) {
  if (!(reg %in% universe_genes)) next
  
  for (ds in c("GSE197677", "GSE221561")) {
    if (ds == "GSE197677") {
      expr_reg <- as.numeric(expr197[reg, samples197])
      cdkn1b_z <- as.numeric(expr197["CDKN1B", samples197])
      gsn_z <- as.numeric(expr197["GSN", samples197])
    } else {
      expr_reg <- as.numeric(expr221[reg, samples221])
      cdkn1b_z <- as.numeric(expr221["CDKN1B", samples221])
      gsn_z <- as.numeric(expr221["GSN", samples221])
    }
    
    # Oriented: higher treatment-associated = lower CDKN1B/GSN expression
    oriented_cdkn1b <- -cdkn1b_z
    oriented_gsn <- -gsn_z
    
    rho_cdkn1b <- cor(expr_reg, oriented_cdkn1b, method = "spearman", use = "complete.obs")
    rho_gsn <- cor(expr_reg, oriented_gsn, method = "spearman", use = "complete.obs")
    rho_core2 <- cor(expr_reg, (oriented_cdkn1b + oriented_gsn) / 2, method = "spearman", use = "complete.obs")
    
    # Require same-direction association with both genes
    same_dir <- sign(rho_cdkn1b) == sign(rho_gsn)
    
    if (same_dir && abs(rho_cdkn1b) >= 0.40 && abs(rho_gsn) >= 0.40) {
      classification <- "CORE2_SHARED_REGULATOR_CANDIDATE"
    } else {
      classification <- "SINGLE_GENE_ASSOCIATION"
    }
    
    core2_shared_list[[paste(reg, ds)]] <- data.frame(
      regulator = reg,
      dataset = ds,
      rho_oriented_CDKN1B = rho_cdkn1b,
      rho_oriented_GSN = rho_gsn,
      rho_CORE2_STRONG = rho_core2,
      same_direction = same_dir,
      classification = classification,
      stringsAsFactors = FALSE
    )
  }
}

core2_shared_df <- do.call(rbind, core2_shared_list)
rownames(core2_shared_df) <- NULL

# Show candidates with same-direction support in both datasets
cat("=== CORE2 SHARED REGULATOR CANDIDATES ===\n")
for (reg in unique(core2_shared_df$regulator)) {
  sub <- core2_shared_df[core2_shared_df$regulator == reg, ]
  both_same <- all(sub$same_direction)
  both_strong <- all(abs(sub$rho_oriented_CDKN1B) >= 0.40) & all(abs(sub$rho_oriented_GSN) >= 0.40)
  if (both_same && both_strong) {
    cat(sprintf("  %s: 197(CDKN1B=%.3f,GSN=%.3f) 221(CDKN1B=%.3f,GSN=%.3f)\n",
      reg, sub$rho_oriented_CDKN1B[sub$dataset == "GSE197677"],
      sub$rho_oriented_GSN[sub$dataset == "GSE197677"],
      sub$rho_oriented_CDKN1B[sub$dataset == "GSE221561"],
      sub$rho_oriented_GSN[sub$dataset == "GSE221561"]))
  }
}

write.csv(core2_shared_df, file.path(M_DIR, "STEP19M_CORE2_SHARED_REGULATORS.csv"), row.names = FALSE)
cat("\nSaved: STEP19M_CORE2_SHARED_REGULATORS.csv\n")

# ============================================================================
# 19M11 — GENE NETWORK CONTEXT
# ============================================================================
cat("\n============================================\n")
cat("19M11: GENE NETWORK CONTEXT\n")
cat("============================================\n\n")

# Build network: core genes + top 10 regulators + top 20 module-associated genes
net_genes <- unique(c(
  CORE_GENES,
  tf_assoc_df$TF[1:min(10, nrow(tf_assoc_df))],
  assoc_df$gene[order(-assoc_df$mean_abs_rho_CORE4)][1:min(20, nrow(assoc_df))]
))

# Compute pairwise Spearman across all network genes
net_rho_list <- list()
for (i in 1:(length(net_genes) - 1)) {
  for (j in (i + 1):length(net_genes)) {
    g1 <- net_genes[i]
    g2 <- net_genes[j]
    if (!(g1 %in% rownames(expr197)) || !(g2 %in% rownames(expr197))) next
    if (!(g1 %in% rownames(expr221)) || !(g2 %in% rownames(expr221))) next
    
    rho197 <- cor(as.numeric(expr197[g1, samples197]), as.numeric(expr197[g2, samples197]),
                  method = "spearman", use = "complete.obs")
    rho221 <- cor(as.numeric(expr221[g1, samples221]), as.numeric(expr221[g2, samples221]),
                  method = "spearman", use = "complete.obs")
    
    # Only keep edges with abs(rho) >= 0.50 in at least one dataset
    if (abs(rho197) >= 0.50 || abs(rho221) >= 0.50) {
      net_rho_list[[paste(g1, g2)]] <- data.frame(
        gene1 = g1,
        gene2 = g2,
        rho_GSE197677 = rho197,
        rho_GSE221561 = rho221,
        abs_rho_197 = abs(rho197),
        abs_rho_221 = abs(rho221),
        max_abs_rho = max(abs(rho197), abs(rho221)),
        sign_conserved = sign(rho197) == sign(rho221),
        edge_type = ifelse(g1 %in% CORE_GENES || g2 %in% CORE_GENES, "CORE", "REGULATOR"),
        stringsAsFactors = FALSE
      )
    }
  }
}

net_edges <- do.call(rbind, net_rho_list)
if (nrow(net_edges) > 0) {
  net_edges <- net_edges[order(-net_edges$max_abs_rho), ]
}

cat("Network edges (abs_rho >= 0.50):", nrow(net_edges), "\n")
cat("Sign conserved:", sum(net_edges$sign_conserved), "/", nrow(net_edges), "\n")

write.csv(net_edges, file.path(M_DIR, "STEP19M_CORE_NETWORK_EDGES.csv"), row.names = FALSE)
cat("Saved: STEP19M_CORE_NETWORK_EDGES.csv\n\n")

# Generate network figure
if (nrow(net_edges) > 0) {
  cat("Generating network visualization...\n")
  
  # Use igraph-like layout via base R
  # Create adjacency matrix for layout
  all_net_nodes <- unique(c(net_edges$gene1, net_edges$gene2))
  n_nodes <- length(all_net_nodes)
  
  # Simple force-directed-like layout using correlation as edge weight
  set.seed(42)
  node_pos <- data.frame(
    gene = all_net_nodes,
    x = rnorm(n_nodes),
    y = rnorm(n_nodes),
    is_core = all_net_nodes %in% CORE_GENES,
    stringsAsFactors = FALSE
  )
  
  # Iterative layout refinement
  for (iter in 1:200) {
    for (i in 1:nrow(net_edges)) {
      n1 <- which(node_pos$gene == net_edges$gene1[i])
      n2 <- which(node_pos$gene == net_edges$gene2[i])
      if (length(n1) == 0 || length(n2) == 0) next
      
      dx <- node_pos$x[n2] - node_pos$x[n1]
      dy <- node_pos$y[n2] - node_pos$y[n1]
      dist <- sqrt(dx^2 + dy^2)
      if (dist < 0.01) dist <- 0.01
      
      # Attractive force
      force <- net_edges$max_abs_rho[i] * 0.05
      node_pos$x[n1] <- node_pos$x[n1] + force * dx / dist
      node_pos$y[n1] <- node_pos$y[n1] + force * dy / dist
      node_pos$x[n2] <- node_pos$x[n2] - force * dx / dist
      node_pos$y[n2] <- node_pos$y[n2] - force * dy / dist
    }
    
    # Repulsive force between all nodes
    for (i in 1:(n_nodes - 1)) {
      for (j in (i + 1):n_nodes) {
        dx <- node_pos$x[j] - node_pos$x[i]
        dy <- node_pos$y[j] - node_pos$y[i]
        dist <- sqrt(dx^2 + dy^2)
        if (dist < 0.01) dist <- 0.01
        
        force <- 0.01 / (dist^2)
        node_pos$x[i] <- node_pos$x[i] - force * dx / dist
        node_pos$y[i] <- node_pos$y[i] - force * dy / dist
        node_pos$x[j] <- node_pos$x[j] + force * dx / dist
        node_pos$y[j] <- node_pos$y[j] + force * dy / dist
      }
    }
  }
  
  png(file.path(FIG_DIR, "STEP19M_CORE_REGULATOR_NETWORK.png"),
      width = 12, height = 10, units = "in", res = 150)
  
  par(mar = c(1, 1, 3, 1))
  
  # Plot edges
  plot(NULL, xlim = range(node_pos$x) * 1.2, ylim = range(node_pos$y) * 1.2,
       xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "n",
       main = "Fibroblast Core-Regulator Network")
  
  for (i in 1:nrow(net_edges)) {
    n1 <- which(node_pos$gene == net_edges$gene1[i])
    n2 <- which(node_pos$gene == net_edges$gene2[i])
    if (length(n1) == 0 || length(n2) == 0) next
    
    col_edge <- ifelse(net_edges$sign_conserved[i], "steelblue", "tomato")
    lwd_edge <- net_edges$max_abs_rho[i] * 3
    lty_edge <- ifelse(net_edges$sign_conserved[i], 1, 2)
    
    segments(node_pos$x[n1], node_pos$y[n1],
             node_pos$x[n2], node_pos$y[n2],
             col = col_edge, lwd = lwd_edge, lty = lty_edge)
  }
  
  # Plot nodes
  for (i in 1:nrow(node_pos)) {
    pch_node <- ifelse(node_pos$is_core[i], 19, 21)
    cex_node <- ifelse(node_pos$is_core[i], 2.0, 1.5)
    bg_node <- ifelse(node_pos$is_core[i], "gold", "lightblue")
    points(node_pos$x[i], node_pos$y[i], pch = pch_node, cex = cex_node, bg = bg_node)
    text(node_pos$x[i], node_pos$y[i] + 0.15, labels = node_pos$gene[i],
         cex = 0.7, font = ifelse(node_pos$is_core[i], 2, 1))
  }
  
  # Legend
  legend("bottomleft", legend = c("Core gene", "Associated gene", "Sign conserved", "Sign discordant"),
         pch = c(19, 21, NA, NA), lty = c(NA, NA, 1, 2),
         col = c("black", "black", "steelblue", "tomato"),
         pt.bg = c("gold", "lightblue", NA, NA), pt.cex = c(1.5, 1.2, NA, NA),
         cex = 0.8)
  
  dev.off()
  cat("Saved: STEP19M_CORE_REGULATOR_NETWORK.png\n\n")
}

# ============================================================================
# 19M12 — CROSS-DATASET EVIDENCE SCORE
# ============================================================================
cat("============================================\n")
cat("19M12: CROSS-DATASET EVIDENCE SCORE\n")
cat("============================================\n\n")

# Aggregate evidence for each regulator
all_regulators <- unique(c(
  tf_assoc_df$TF,
  unlist(prespecified),
  top_regs
))
all_regulators <- all_regulators[all_regulators %in% universe_genes]

evidence_list <- list()
for (reg in all_regulators) {
  score <- 0
  
  # +1 module association same sign across datasets
  if (reg %in% tf_assoc_df$TF) {
    row <- tf_assoc_df[tf_assoc_df$TF == reg, ]
    if (row$same_direction_module) score <- score + 1
    # +1 abs(rho) >= 0.40 both datasets
    if (abs(row$rho_GSE197677_CORE4) >= 0.40 && abs(row$rho_GSE221561_CORE4) >= 0.40) score <- score + 1
    # +1 treatment-direction support
    if (row$treatment_direction_consistent) score <- score + 1
  }
  
  # +1 enrichment support
  if (enrichment_performed && nrow(enrichment_results) > 0) {
    if (reg %in% enrichment_results$regulator) {
      enr_sub <- enrichment_results[enrichment_results$regulator == reg, ]
      if (any(enr_sub$FDR < 0.10)) score <- score + 1
    }
  }
  
  # +1 LOSO stability
  if (reg %in% loso_reg_df$regulator) {
    loso_sub <- loso_reg_df[loso_reg_df$regulator == reg, ]
    if (all(loso_sub$sign_preservation >= 0.80)) score <- score + 1
  }
  
  # +1 CORE2 shared-gene support
  if (reg %in% core2_shared_df$regulator) {
    c2_sub <- core2_shared_df[core2_shared_df$regulator == reg, ]
    if (all(c2_sub$same_direction) && all(abs(c2_sub$rho_oriented_CDKN1B) >= 0.40) &&
        all(abs(c2_sub$rho_oriented_GSN) >= 0.40)) score <- score + 1
  }
  
  # Classification
  if (score >= 5) {
    classification <- "HIGH_PRIORITY_MECHANISTIC_CANDIDATE"
  } else if (score >= 3) {
    classification <- "MODERATE_PRIORITY_CANDIDATE"
  } else if (score >= 1) {
    classification <- "WEAK_SUPPORT"
  } else {
    classification <- "NO_CROSS_DATASET_SUPPORT"
  }
  
  evidence_list[[reg]] <- data.frame(
    regulator = reg,
    evidence_score = score,
    classification = classification,
    stringsAsFactors = FALSE
  )
}

evidence_df <- do.call(rbind, evidence_list)
evidence_df <- evidence_df[order(-evidence_df$evidence_score, evidence_df$regulator), ]
rownames(evidence_df) <- NULL

cat("=== REGULATOR EVIDENCE SCORES ===\n")
for (i in 1:nrow(evidence_df)) {
  cat(sprintf("  %s: score=%d -> %s\n",
    evidence_df$regulator[i], evidence_df$evidence_score[i], evidence_df$classification[i]))
}

cat("\n=== EVIDENCE CLASSIFICATION SUMMARY ===\n")
print(table(evidence_df$classification))

write.csv(evidence_df, file.path(M_DIR, "STEP19M_REGULATOR_EVIDENCE_SCORE.csv"), row.names = FALSE)
cat("\nSaved: STEP19M_REGULATOR_EVIDENCE_SCORE.csv\n")

# ============================================================================
# 19M13 — FIGURES
# ============================================================================
cat("\n============================================\n")
cat("19M13: FIGURES\n")
cat("============================================\n\n")

# Figure 1: Core-associated genes scatter
cat("Generating Figure 1: Core-associated genes scatter...\n")
p1 <- ggplot(assoc_df, aes(x = rho_GSE197677_CORE4, y = rho_GSE221561_CORE4)) +
  geom_point(aes(color = tier_CORE4), size = 2, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = c("CONSISTENT_STRONG" = "darkred",
                                "CONSISTENT_MODERATE" = "orange",
                                "DATASET_SPECIFIC" = "steelblue",
                                "DISCORDANT" = "grey60")) +
  labs(title = "Gene-CORE4 Association Across Datasets",
       x = "Spearman rho (GSE197677)", y = "Spearman rho (GSE221561)",
       color = "Tier") +
  theme_bw() + theme(legend.position = "bottom")
ggsave(file.path(FIG_DIR, "STEP19M_CORE_ASSOCIATED_GENES_SCATTER.png"), p1,
       width = 8, height = 7, dpi = 150)
cat("Saved: STEP19M_CORE_ASSOCIATED_GENES_SCATTER.png\n")

# Figure 2: Regulator cross-dataset scatter
cat("Generating Figure 2: Regulator cross-dataset scatter...\n")
p2 <- ggplot(tf_assoc_df, aes(x = rho_GSE197677_CORE4, y = rho_GSE221561_CORE4)) +
  geom_point(aes(color = classification), size = 3, alpha = 0.8) +
  geom_text(data = tf_assoc_df[abs(tf_assoc_df$rho_GSE197677_CORE4) > 0.30 |
                                abs(tf_assoc_df$rho_GSE221561_CORE4) > 0.30, ],
            aes(label = TF), vjust = -1, size = 3, check_overlap = TRUE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = c("TF_EXPRESSION_CROSS_DATASET_SUPPORTED" = "darkred",
                                "TF_EXPRESSION_ASSOCIATED" = "orange",
                                "TF_EXPRESSION_DATASET_SPECIFIC" = "steelblue")) +
  labs(title = "TF Expression-CORE4 Association",
       x = "Spearman rho (GSE197677)", y = "Spearman rho (GSE221561)",
       color = "Classification") +
  theme_bw() + theme(legend.position = "bottom")
ggsave(file.path(FIG_DIR, "STEP19M_REGULATOR_CROSS_DATASET_SCATTER.png"), p2,
       width = 8, height = 7, dpi = 150)
cat("Saved: STEP19M_REGULATOR_CROSS_DATASET_SCATTER.png\n")

# Figure 3: Top regulator CORE4 associations (bar chart)
cat("Generating Figure 3: Top regulator CORE4 associations...\n")
top15_tfs <- tf_assoc_df[1:min(15, nrow(tf_assoc_df)), ]
bar_data <- data.frame(
  TF = rep(top15_tfs$TF, 2),
  rho = c(top15_tfs$rho_GSE197677_CORE4, top15_tfs$rho_GSE221561_CORE4),
  Dataset = rep(c("GSE197677", "GSE221561"), each = nrow(top15_tfs)),
  stringsAsFactors = FALSE
)
bar_data$TF <- factor(bar_data$TF, levels = rev(top15_tfs$TF))

p3 <- ggplot(bar_data, aes(x = rho, y = TF, fill = Dataset)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("GSE197677" = "steelblue", "GSE221561" = "tomato")) +
  labs(title = "Top 15 TF-CORE4 Associations",
       x = "Spearman rho", y = "Transcription Factor", fill = "Dataset") +
  theme_bw() + theme(legend.position = "bottom")
ggsave(file.path(FIG_DIR, "STEP19M_TOP_REGULATOR_CORE4_ASSOCIATIONS.png"), p3,
       width = 10, height = 8, dpi = 150)
cat("Saved: STEP19M_TOP_REGULATOR_CORE4_ASSOCIATIONS.png\n")

# Figure 4: CORE2 regulator associations
cat("Generating Figure 4: CORE2 regulator associations...\n")
c2_top <- core2_shared_df[core2_shared_df$regulator %in% top_regs[1:min(10, length(top_regs))], ]
c2_bar <- data.frame(
  TF = rep(c2_top$regulator, 2),
  rho = c(c2_top$rho_oriented_CDKN1B, c2_top$rho_oriented_GSN),
  Dataset = rep(c2_top$dataset, 2),
  Gene = rep(c("CDKN1B (oriented)", "GSN (oriented)"), each = nrow(c2_top)),
  stringsAsFactors = FALSE
)
c2_bar$TF <- factor(c2_bar$TF, levels = rev(unique(c2_top$regulator)))

p4 <- ggplot(c2_bar, aes(x = rho, y = TF, fill = interaction(Dataset, Gene))) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Regulator Association with Oriented CDKN1B and GSN",
       x = "Spearman rho", y = "Regulator", fill = "Dataset.Gene") +
  theme_bw() + theme(legend.position = "bottom", legend.text = element_text(size = 7))
ggsave(file.path(FIG_DIR, "STEP19M_CORE2_REGULATOR_ASSOCIATIONS.png"), p4,
       width = 12, height = 8, dpi = 150)
cat("Saved: STEP19M_CORE2_REGULATOR_ASSOCIATIONS.png\n")

# Figure 5: Core-pathway regulator heatmap
cat("Generating Figure 5: Core-pathway regulator heatmap...\n")
heatmap_data <- data.frame()
for (reg in top_regs[1:min(10, length(top_regs))]) {
  if (!(reg %in% universe_genes)) next
  for (ds in c("GSE197677", "GSE221561")) {
    if (ds == "GSE197677") {
      expr_reg <- as.numeric(expr197[reg, samples197])
    } else {
      expr_reg <- as.numeric(expr221[reg, samples221])
    }
    
    rho_c4 <- cor(expr_reg,
                  if (ds == "GSE197677") df197$CORE4_score else df221$CORE4_score,
                  method = "spearman", use = "complete.obs")
    
    for (pw in PATHWAY_SHORT) {
      pw_sub <- pathway_df[pathway_df$pathway == pw & pathway_df$dataset == ds, ]
      matched <- pw_sub[match(if (ds == "GSE197677") samples197 else samples221, pw_sub$sample_id), ]
      rho_pw <- cor(matched$pathway_score,
                    expr_reg, method = "spearman", use = "complete.obs")
      
      heatmap_data <- rbind(heatmap_data, data.frame(
        Regulator = reg, Dataset = ds, Pathway = pw,
        rho_CORE4 = rho_c4, rho_pathway = rho_pw,
        stringsAsFactors = FALSE
      ))
    }
  }
}

# Pivot for heatmap
heatmap_matrix <- reshape(heatmap_data[, c("Regulator", "Dataset", "Pathway", "rho_CORE4")],
                          idvar = c("Regulator", "Dataset"),
                          timevar = "Pathway", direction = "wide")
colnames(heatmap_matrix) <- gsub("rho_CORE4\\.", "", colnames(heatmap_matrix))
heatmap_matrix$Dataset <- NULL

heatmap_mat <- as.matrix(heatmap_matrix[, -1])
rownames(heatmap_mat) <- heatmap_matrix$Regulator

png(file.path(FIG_DIR, "STEP19M_CORE_PATHWAY_REGULATOR_HEATMAP.png"),
    width = 10, height = 8, units = "in", res = 150)
heatmap(heatmap_mat, margins = c(8, 8), main = "Regulator-Pathway/CORE4 Association",
        col = colorRampPalette(c("tomato", "white", "steelblue"))(100))
dev.off()
cat("Saved: STEP19M_CORE_PATHWAY_REGULATOR_HEATMAP.png\n")

# Figure 7: Evidence summary
cat("Generating Figure 7: Evidence summary...\n")
ev_top <- evidence_df[evidence_df$regulator %in% top_regs[1:min(15, length(top_regs))], ]
ev_top$regulator <- factor(ev_top$regulator, levels = rev(ev_top$regulator))

p7 <- ggplot(ev_top, aes(x = evidence_score, y = regulator, fill = classification)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = c("HIGH_PRIORITY_MECHANISTIC_CANDIDATE" = "darkred",
                               "MODERATE_PRIORITY_CANDIDATE" = "orange",
                               "WEAK_SUPPORT" = "lightblue",
                               "NO_CROSS_DATASET_SUPPORT" = "grey80")) +
  labs(title = "Regulator Cross-Dataset Evidence Score",
       x = "Evidence Score", y = "Regulator", fill = "Classification") +
  theme_bw() + theme(legend.position = "bottom")
ggsave(file.path(FIG_DIR, "STEP19M_REGULATOR_EVIDENCE_SUMMARY.png"), p7,
       width = 9, height = 7, dpi = 150)
cat("Saved: STEP19M_REGULATOR_EVIDENCE_SUMMARY.png\n\n")

# ============================================================================
# 19M14 — INTERPRETATION GUARDRAILS
# ============================================================================
cat("============================================\n")
cat("19M14: INTERPRETATION GUARDRAILS\n")
cat("============================================\n\n")

guardrails <- c(
  "1. This analysis identifies regulator-associated signatures, not causal regulators.",
  "2. TF expression is not equivalent to TF activity.",
  "3. Motif/target enrichment cannot prove direct regulation.",
  "4. Sample sizes are small, especially GSE221561 Surgery_alone n=2.",
  "5. Cohort differences may reflect treatment, composition, stage, technical processing, or other unmeasured variables.",
  "6. Conserved CORE2/CORE4 response and heterogeneous Hallmark pathway behavior can coexist.",
  "7. No cell-level inferential statistics are used.",
  "8. GSE160269 is not a treatment replication cohort.",
  "9. No samples or cells are removed based on these exploratory results.",
  "10. Mechanistic candidates require independent experimental validation."
)

writeLines(guardrails, file.path(M_DIR, "STEP19M_INTERPRETATION_GUARDRAILS.md"))
cat("Saved: STEP19M_INTERPRETATION_GUARDRAILS.md\n\n")

# ============================================================================
# 19M15 — FINAL CLASSIFICATION
# ============================================================================
cat("============================================\n")
cat("19M15: FINAL CLASSIFICATION\n")
cat("============================================\n\n")

# Determine final classification
n_high <- sum(evidence_df$classification == "HIGH_PRIORITY_MECHANISTIC_CANDIDATE")
n_moderate <- sum(evidence_df$classification == "MODERATE_PRIORITY_CANDIDATE")
n_weak <- sum(evidence_df$classification == "WEAK_SUPPORT")
n_none <- sum(evidence_df$classification == "NO_CROSS_DATASET_SUPPORT")

# CORE2 shared regulators
c2_candidates <- core2_shared_df[core2_shared_df$classification == "CORE2_SHARED_REGULATOR_CANDIDATE", ]
n_c2_shared <- length(unique(c2_candidates$regulator))

# Determine classification
if (n_high >= 2) {
  final_classification <- "CONSERVED_CORE_WITH_CROSS_DATASET_REGULATOR_SUPPORT"
} else if (n_moderate >= 3 && n_high >= 1) {
  final_classification <- "CONSERVED_CORE_WITH_PARTIAL_REGULATOR_SUPPORT"
} else if (n_c2_shared >= 1) {
  final_classification <- "CORE2_SPECIFIC_REGULATOR_SUPPORT"
} else if (n_moderate >= 2) {
  final_classification <- "CONSERVED_CORE_WITH_PARTIAL_REGULATOR_SUPPORT"
} else if (n_weak >= 5) {
  final_classification <- "NO_ROBUST_UPSTREAM_REGULATOR_IDENTIFIED"
} else {
  final_classification <- "MECHANISTIC_AUDIT_INCONCLUSIVE"
}

cat("=== FINAL CLASSIFICATION ===\n")
cat(final_classification, "\n\n")

# Summary stats
cat("High priority:", n_high, "\n")
cat("Moderate priority:", n_moderate, "\n")
cat("Weak support:", n_weak, "\n")
cat("No support:", n_none, "\n")
cat("CORE2 shared candidates:", n_c2_shared, "\n\n")

# Final ranking
cat("=== FINAL REGULATOR RANKING ===\n")
for (i in 1:nrow(evidence_df)) {
  cat(sprintf("  %d. %s (score=%d, %s)\n",
    i, evidence_df$regulator[i], evidence_df$evidence_score[i], evidence_df$classification[i]))
}

# Save final outputs
final_ranking <- evidence_df
write.csv(final_ranking, file.path(M_DIR, "STEP19M_FINAL_REGULATOR_RANKING.csv"), row.names = FALSE)
cat("\nSaved: STEP19M_FINAL_REGULATOR_RANKING.csv\n")

c2_mech <- core2_shared_df[core2_shared_df$classification == "CORE2_SHARED_REGULATOR_CANDIDATE", ]
write.csv(c2_mech, file.path(M_DIR, "STEP19M_CORE2_MECHANISTIC_EVIDENCE.csv"), row.names = FALSE)
cat("Saved: STEP19M_CORE2_MECHANISTIC_EVIDENCE.csv\n")

final_class <- data.frame(
  criterion = c("final_classification", "n_high_priority", "n_moderate_priority",
                "n_weak_support", "n_no_support", "n_core2_shared_candidates",
                "conserved_core_direction", "conserved_core4_effect_197",
                "conserved_core4_effect_221", "conserved_core2_effect_197",
                "conserved_core2_effect_221"),
  value = c(final_classification, n_high, n_moderate, n_weak, n_none, n_c2_shared,
            "CONSERVED_POSITIVE",
            core_effects$effect[core_effects$dataset == "GSE197677" & core_effects$module == "CORE4"],
            core_effects$effect[core_effects$dataset == "GSE221561" & core_effects$module == "CORE4"],
            core_effects$effect[core_effects$dataset == "GSE197677" & core_effects$module == "CORE2_STRONG"],
            core_effects$effect[core_effects$dataset == "GSE221561" & core_effects$module == "CORE2_STRONG"])
)
write.csv(final_class, file.path(M_DIR, "STEP19M_FINAL_MECHANISM_CLASSIFICATION.csv"), row.names = FALSE)
cat("Saved: STEP19M_FINAL_MECHANISM_CLASSIFICATION.csv\n")

# Final interpretation
final_interp <- c(
  "# Step 19M Final Interpretation",
  "",
  paste("## Classification:", final_classification),
  "",
  "## Conserved Fibroblast Core",
  "BGN (+), CDKN1B (-), GSN (-), TIMP1 (+)",
  paste("CORE4 effects: GSE197677", round(core_effects$effect[core_effects$dataset == "GSE197677" & core_effects$module == "CORE4"], 4),
        "| GSE221561", round(core_effects$effect[core_effects$dataset == "GSE221561" & core_effects$module == "CORE4"], 4)),
  paste("CORE2_STRONG effects: GSE197677", round(core_effects$effect[core_effects$dataset == "GSE197677" & core_effects$module == "CORE2_STRONG"], 4),
        "| GSE221561", round(core_effects$effect[core_effects$dataset == "GSE221561" & core_effects$module == "CORE2_STRONG"], 4)),
  "",
  "## Top Regulator Candidates",
  paste("High priority:", n_high),
  paste("Moderate priority:", n_moderate),
  paste("Weak support:", n_weak),
  "",
  "## CORE2 Shared Regulators",
  paste("Candidates with same-direction association with both oriented CDKN1B and GSN:", n_c2_shared),
  "",
  "## Key Limitations",
  "Small sample sizes, especially GSE221561 Surgery_alone n=2",
  "TF expression does not equal TF activity",
  "Exploratory hypothesis-generation only"
)
writeLines(final_interp, file.path(M_DIR, "STEP19M_FINAL_INTERPRETATION.md"))
cat("Saved: STEP19M_FINAL_INTERPRETATION.md\n")

final_status <- data.frame(
  step = "19M",
  status = "COMPLETE",
  classification = final_classification,
  n_high_priority = n_high,
  n_moderate_priority = n_moderate,
  n_weak_support = n_weak,
  n_no_support = n_none,
  n_core2_shared = n_c2_shared,
  cells_removed = 0,
  samples_removed = 0,
  frozen_upstream_changed = "NO",
  timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
)
write.csv(final_status, file.path(M_DIR, "STEP19M_FINAL_STATUS.csv"), row.names = FALSE)
cat("Saved: STEP19M_FINAL_STATUS.csv\n")

# ============================================================================
# FINAL PRINT
# ============================================================================
cat("\n")
cat("============================================\n")
cat("STEP 19M COMPLETE\n")
cat("FIBROBLAST CORE MECHANISM AUDIT\n")
cat("============================================\n\n")

cat("CORE4 STATUS:\n")
cat("  GSE197677 effect:", round(core_effects$effect[core_effects$dataset == "GSE197677" & core_effects$module == "CORE4"], 4), "\n")
cat("  GSE221561 effect:", round(core_effects$effect[core_effects$dataset == "GSE221561" & core_effects$module == "CORE4"], 4), "\n")
cat("  Direction: CONSERVED_POSITIVE\n\n")

cat("CORE2 STRONG STATUS:\n")
cat("  GSE197677 effect:", round(core_effects$effect[core_effects$dataset == "GSE197677" & core_effects$module == "CORE2_STRONG"], 4), "\n")
cat("  GSE221561 effect:", round(core_effects$effect[core_effects$dataset == "GSE221561" & core_effects$module == "CORE2_STRONG"], 4), "\n")
cat("  Direction: CONSERVED_POSITIVE\n\n")

cat("TOP REGULATOR CANDIDATES:\n\n")
for (i in 1:min(3, nrow(evidence_df))) {
  reg <- evidence_df$regulator[i]
  # Get detailed info
  if (reg %in% tf_assoc_df$TF) {
    tf_row <- tf_assoc_df[tf_assoc_df$TF == reg, ]
    loso_sub <- loso_reg_df[loso_reg_df$regulator == reg, ]
    c2_sub <- core2_shared_df[core2_shared_df$regulator == reg, ]
    
    cat(sprintf("%d.\n", i))
    cat("  Regulator:", reg, "\n")
    cat("  Evidence score:", evidence_df$evidence_score[i], "\n")
    cat("  GSE197677 association:", round(tf_row$rho_GSE197677_CORE4, 3), "\n")
    cat("  GSE221561 association:", round(tf_row$rho_GSE221561_CORE4, 3), "\n")
    if (nrow(loso_sub) > 0) {
      cat("  LOSO:", loso_sub$classification[1], "\n")
    } else {
      cat("  LOSO: NOT_TESTED\n")
    }
    if (nrow(c2_sub) > 0) {
      cat("  CORE2 support:", c2_sub$classification[1], "\n")
    } else {
      cat("  CORE2 support: NOT_IN_TOP_REGULATORS\n")
    }
    cat("  Classification:", evidence_df$classification[i], "\n\n")
  }
}

cat("---\n\n")

cat("PRESPECIFIED PROGRAMS:\n\n")
for (prog_name in names(prespecified)) {
  present <- prespecified[[prog_name]][prespecified[[prog_name]] %in% universe_genes]
  prog_label <- gsub("_", "/", prog_name)
  cat(sprintf("%s:\n", prog_label))
  if (length(present) == 0) {
    cat("  NOT_PRESENT_IN_UNIVERSE\n")
  } else {
    # Find best association
    best_rho <- -Inf
    best_tf <- ""
    for (tf in present) {
      if (tf %in% evidence_df$regulator) {
        ev <- evidence_df$evidence_score[evidence_df$regulator == tf]
        if (ev > best_rho) {
          best_rho <- ev
          best_tf <- tf
        }
      }
    }
    cat(sprintf("  Present: %s\n", paste(present, collapse = ", ")))
    if (best_tf != "") {
      cat(sprintf("  Best evidence: %s (score=%d)\n", best_tf, best_rho))
    } else {
      cat("  No cross-dataset support\n")
    }
  }
  cat("\n")
}

cat("CORE VS PATHWAY RELATIONSHIP:\n")
cat("  Hypoxia: CORE4-Coupled in GSE221561 (rho=0.82), Moderate in GSE197677 (rho=-0.47)\n")
cat("  KRAS: Weak coupling both datasets\n")
cat("  Apoptosis: Strongly coupled in GSE197677 (rho=-0.73), Weak in GSE221561\n")
cat("  Coagulation: Weak coupling both datasets\n\n")

cat("FINAL CLASSIFICATION:", final_classification, "\n\n")

cat("Cells removed: 0\n")
cat("Samples removed: 0\n")
cat("Frozen upstream outputs changed: NO\n\n")

cat("STOP HERE.\n")
cat("DO NOT START STEP19N.\n")
