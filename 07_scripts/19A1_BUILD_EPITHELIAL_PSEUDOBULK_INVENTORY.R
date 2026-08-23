#!/usr/bin/env Rscript
# ============================================================
# STEP 19A1: BUILD CANONICAL TUMOR EPITHELIAL SAMPLE INVENTORY
# ============================================================

library(data.table)

cat("============================================================\n")
cat("STEP 19A1: TUMOR EPITHELIAL SAMPLE INVENTORY\n")
cat("============================================================\n\n")

OUTDIR <- "06_tables/cross_dataset/epithelial_pseudobulk"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 1. READ STEP 18F CAPACITY
# ============================================================

cat("1. Reading Step 18F capacity...\n\n")

capacity <- fread("06_tables/cross_dataset/epithelial_analysis/STEP18F_three_dataset_epithelial_analysis_capacity.csv")
print(capacity)

# ============================================================
# 2. READ STEP 18F GSE160269 INVENTORY
# ============================================================

cat("\n2. Reading Step 18F GSE160269 inventory...\n\n")

g160269_inv <- fread("06_tables/cross_dataset/epithelial_analysis/STEP18F_GSE160269_tumor_epithelial_sample_inventory.csv")
cat("GSE160269 tumor samples:", nrow(g160269_inv), "\n")

# ============================================================
# 3. READ PRIMARY LINEAGE COUNTS
# ============================================================

cat("\n3. Reading primary lineage counts...\n\n")

counts <- fread("06_tables/cross_dataset/STEP16_sample_primary_lineage_cell_counts_WIDE.csv")

# ============================================================
# 4. BUILD GSE221561 INVENTORY
# ============================================================

cat("\n4. Building GSE221561 inventory...\n\n")

g221 <- counts[dataset == "GSE221561" & tissue == "Tumor"]
cat("GSE221561 tumor samples:", nrow(g221), "\n")

g221_inv <- data.table(
  dataset = "GSE221561",
  sample = g221$sample_id,
  patient = g221$patient_id,
  tissue = "Tumor",
  treatment_original = g221$treatment_group,
  treatment_standardized = g221$treatment_group,
  epithelial_cells = g221$Epithelial,
  raw_counts_available = NA,
  raw_gene_features = NA,
  gene_id_format = "TO_BE_CHECKED",
  eligible_ge30 = g221$Epithelial >= 30,
  eligible_ge100 = g221$Epithelial >= 100,
  analysis_role = "PRIMARY_TREATMENT_RESPONSE",
  internal_contrast_group = fifelse(
    g221$treatment_group == "Neoadjuvant_treated",
    "Neoadjuvant_treated",
    "Surgery_alone"
  )
)

# Mark small reference group
g221_inv[treatment_standardized == "Surgery_alone",
         internal_contrast_group := paste0(internal_contrast_group, " [SMALL_REFERENCE_GROUP n=2]")]
g221_inv[, internal_contrast_group := gsub(" \\[SMALL_REFERENCE_GROUP n=2\\]$",
  ifelse(treatment_standardized == "Surgery_alone", " [SMALL_REFERENCE_GROUP n=2]", ""),
  internal_contrast_group)]

# Fix: just set directly
g221_inv[treatment_standardized == "Surgery_alone",
         internal_contrast_group := "Surgery_alone [SMALL_REFERENCE_GROUP n=2]"]
g221_inv[treatment_standardized == "Neoadjuvant_treated",
         internal_contrast_group := "Neoadjuvant_treated"]

cat("GSE221561 groups:\n")
print(table(g221_inv$treatment_standardized))

# ============================================================
# 5. BUILD GSE197677 INVENTORY
# ============================================================

cat("\n5. Building GSE197677 inventory...\n\n")

g197_broad <- fread("06_tables/cross_dataset/STEP16A_GSE197677_sample_by_broad_celltype_REBUILT.csv")
g197_tumor_epi <- g197_broad[tissue == "Tumor" & broad_celltype == "Epithelial"]
g197_epith_per_sample <- g197_tumor_epi[, .(epithelial_cells = sum(n_cells)), by = .(patient_id, sample_id, treatment_group)]

cat("GSE197677 tumor samples with epithelial:", nrow(g197_epith_per_sample), "\n")

g197_inv <- data.table(
  dataset = "GSE197677",
  sample = g197_epith_per_sample$sample_id,
  patient = g197_epith_per_sample$patient_id,
  tissue = "Tumor",
  treatment_original = g197_epith_per_sample$treatment_group,
  treatment_standardized = g197_epith_per_sample$treatment_group,
  epithelial_cells = g197_epith_per_sample$epithelial_cells,
  raw_counts_available = NA,
  raw_gene_features = NA,
  gene_id_format = "TO_BE_CHECKED",
  eligible_ge30 = g197_epith_per_sample$epithelial_cells >= 30,
  eligible_ge100 = g197_epith_per_sample$epithelial_cells >= 100,
  analysis_role = "PRIMARY_INTERNAL_REPLICATION",
  internal_contrast_group = fifelse(
    g197_epith_per_sample$treatment_group == "Neoadjuvant_chemotherapy",
    "NACT",
    "nNACT"
  )
)

cat("GSE197677 groups:\n")
print(table(g197_inv$internal_contrast_group))

# ============================================================
# 6. BUILD GSE160269 INVENTORY
# ============================================================

cat("\n6. Building GSE160269 inventory...\n\n")

g160269_inv_full <- data.table(
  dataset = "GSE160269",
  sample = g160269_inv$sample,
  patient = g160269_inv$patient,
  tissue = "Tumor",
  treatment_original = "Untreated_baseline",
  treatment_standardized = "Untreated_baseline",
  epithelial_cells = g160269_inv$epithelial_cells,
  raw_counts_available = NA,
  raw_gene_features = NA,
  gene_id_format = "TO_BE_CHECKED",
  eligible_ge30 = g160269_inv$primary_epithelial_analysis_eligible,
  eligible_ge100 = g160269_inv$epithelial_cells >= 100,
  analysis_role = "UNTREATED_BASELINE_REFERENCE",
  internal_contrast_group = "NONE"
)

cat("GSE160269 eligible_ge30:", sum(g160269_inv_full$eligible_ge30), "\n")
cat("GSE160269 eligible_ge100:", sum(g160269_inv_full$eligible_ge100), "\n")

# ============================================================
# 7. COMBINE AND SAVE
# ============================================================

cat("\n7. Combining inventory...\n\n")

inventory <- rbind(g221_inv, g197_inv, g160269_inv_full)

fwrite(inventory, file.path(OUTDIR, "STEP19A_tumor_epithelial_sample_inventory.csv"))

cat("Total samples:", nrow(inventory), "\n")
cat("  GSE221561:", sum(inventory$dataset == "GSE221561"), "\n")
cat("  GSE197677:", sum(inventory$dataset == "GSE197677"), "\n")
cat("  GSE160269:", sum(inventory$dataset == "GSE160269"), "\n")
cat("\nEligible (ge30):", sum(inventory$eligible_ge30), "\n")
cat("Eligible (ge100):", sum(inventory$eligible_ge100), "\n")
cat("\nSaved: STEP19A_tumor_epithelial_sample_inventory.csv\n")
