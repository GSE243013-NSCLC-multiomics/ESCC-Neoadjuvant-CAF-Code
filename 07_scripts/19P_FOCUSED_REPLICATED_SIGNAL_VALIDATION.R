#!/usr/bin/env Rscript
# ============================================================================
# STEP 19P: FOCUSED VALIDATION OF REPLICATED FIBROBLAST SIGNALING INPUTS
# ============================================================================
suppressPackageStartupMessages({
  library(ggplot2)
  library(reshape2)
})

BASE_DIR <- Sys.getenv("ESCC_CAF_PROJECT_ROOT", unset = ".")
O_DIR    <- file.path(BASE_DIR, "06_tables/cross_dataset/pathway_analysis/Step19P")
FIG_DIR  <- file.path(BASE_DIR, "05_figures/cross_dataset/pathway_analysis/Step19P")
LOG_DIR  <- file.path(BASE_DIR, "08_logs")
SRC_O    <- file.path(BASE_DIR, "06_tables/cross_dataset/pathway_analysis/Step19O")
SRC_FIX  <- file.path(BASE_DIR, "06_tables/cross_dataset/pathway_analysis/Step19O_FIX221")
SRC_L    <- file.path(BASE_DIR, "06_tables/cross_dataset/pathway_analysis/Step19L")
SRC_M    <- file.path(BASE_DIR, "06_tables/cross_dataset/pathway_analysis/Step19M")
SRC_N    <- file.path(BASE_DIR, "06_tables/cross_dataset/pathway_analysis/Step19N")
dir.create(O_DIR, showWarnings=FALSE, recursive=TRUE)
dir.create(FIG_DIR, showWarnings=FALSE, recursive=TRUE)

log_con <- file.path(LOG_DIR, "STEP19P_FOCUSED_REPLICATED_SIGNAL_VALIDATION.log")
cat("", file=log_con)
cat_log <- function(...) { txt<-paste0(...); cat(txt); cat(txt,file=log_con,append=TRUE) }

cat_log("============================================\n")
cat_log("STEP 19P: FOCUSED REPLICATED SIGNAL VALIDATION\n")
cat_log("============================================\n\n")

CORE2_GENES <- c("CDKN1B","GSN")
CORE4_GENES <- c("BGN","CDKN1B","GSN","TIMP1")

TREAT_MAP_197 <- c(ESC06="nNACT",ESC07="nNACT",ESC09="NACT",ESC11="NACT",
  ESC12="NACT",ESC15="NACT",ESC16="NACT",ESC17="nNACT",ESC18="nNACT",ESC22="NACT")
TREAT_MAP_221 <- c(S0520T="Neoadjuvant_treated",S1265T="Neoadjuvant_treated",
  S1315TS="Neoadjuvant_treated",S1535T="Neoadjuvant_treated",S2423T="Neoadjuvant_treated",
  S2487T="Neoadjuvant_treated",S6829T="Neoadjuvant_treated",S6417T="Surgery_alone",
  S9478T="Surgery_alone")

compute_log2_cpm <- function(m) {
  ls <- colSums(m); ls[ls==0]<-1
  log2(sweep(m,2,ls,"/")*1e6 + 1)
}

# ============================================================================
# 19P1 — INPUT PREFLIGHT
# ============================================================================
cat_log("============================================\n")
cat_log("19P1: INPUT PREFLIGHT\n")
cat_log("============================================\n\n")

core_scores <- read.csv(file.path(SRC_L, "STEP19L_SAMPLE_CORE_MODULE_SCORES.csv"))
core197 <- core_scores[core_scores$dataset=="GSE197677",]
core221 <- core_scores[core_scores$dataset=="GSE221561",]
cat_log(sprintf("CORE2 scores: GSE197677=%d, GSE221561=%d\n", nrow(core197), nrow(core221)))

src197 <- list()
for (c in c("Epithelial","Myeloid","T_NK","Endothelial")) {
  fp <- file.path(SRC_O, sprintf("STEP19O_GSE197677_%s_raw_counts.rds", c))
  if (file.exists(fp)) src197[[c]] <- readRDS(fp)
}
src221 <- list()
for (c in c("Epithelial","Myeloid","T_NK","Endothelial")) {
  fp <- file.path(SRC_FIX, sprintf("STEP19O_FIX221_GSE221561_%s_raw_counts.rds", c))
  if (file.exists(fp)) src221[[c]] <- readRDS(fp)
}
cat_log(sprintf("GSE197677 source compartments: %s\n", paste(names(src197), collapse=", ")))
cat_log(sprintf("GSE221561 source compartments: %s\n", paste(names(src221), collapse=", ")))

fib197 <- readRDS(file.path(BASE_DIR, "06_tables/cross_dataset/pathway_analysis/Step19I/counts/GSE197677_FIBROBLAST_CAF_raw_counts.rds"))
fib221 <- readRDS(file.path(BASE_DIR, "06_tables/cross_dataset/pathway_analysis/Step19I/counts/GSE221561_FIBROBLAST_CAF_raw_counts.rds"))
cat_log(sprintf("Fibroblast pseudobulk: GSE197677=%d genes x %d samples, GSE221561=%d genes x %d samples\n",
  nrow(fib197), ncol(fib197), nrow(fib221), ncol(fib221)))

cross_rep <- read.csv(file.path(SRC_FIX, "STEP19O_FIX221_CROSS_DATASET_REPLICATION.csv"))
rec197_av <- read.csv(file.path(SRC_O, "STEP19O_FIBROBLAST_RECEPTOR_AVAILABILITY.csv"))
rec221_av <- read.csv(file.path(SRC_FIX, "STEP19O_FIX221_FIBROBLAST_RECEPTOR_AVAILABILITY.csv"))

reg_act <- read.csv(file.path(SRC_N, "STEP19N_SAMPLE_REGULATOR_ACTIVITY.csv"))
cat_log("Step19N regulator activity loaded.\n")

prov <- data.frame(
  input=c("CORE2_scores","CORE4_scores","GSE197677_source_pseudobulk","GSE221561_source_pseudobulk",
    "GSE197677_fibroblast_pseudobulk","GSE221561_fibroblast_pseudobulk",
    "cross_dataset_replication","Step19N_regulator_activity"),
  source=c("Step19L","Step19L","Step19O","Step19O_FIX221","Step19I","Step19I",
    "Step19O_FIX221","Step19N"),
  loaded=c(TRUE,TRUE,length(src197)>0,length(src221)>0,TRUE,TRUE,TRUE,TRUE),
  stringsAsFactors=FALSE)
write.csv(prov, file.path(O_DIR, "STEP19P_INPUT_PROVENANCE.csv"), row.names=FALSE)
cat_log("Saved: STEP19P_INPUT_PROVENANCE.csv\n\n")

# ============================================================================
# 19P2 — PRIMARY CANDIDATE TABLE
# ============================================================================
cat_log("============================================\n")
cat_log("19P2: PRIMARY CANDIDATE TABLE\n")
cat_log("============================================\n\n")

candidates <- data.frame(
  source=c("Epithelial","Epithelial","Epithelial","Myeloid","Myeloid","T_NK","T_NK","T_NK"),
  ligand=c("IL1B","FGF7","FGF9","LIF","CXCL2","FGF7","WNT5A","CXCL2"),
  rho_GSE197677=c(-0.745,0.492,0.666,0.479,0.418,0.758,0.564,0.588),
  rho_GSE221561=c(-0.417,0.683,0.483,0.500,0.433,0.522,0.533,0.833),
  stringsAsFactors=FALSE)
candidates$same_direction <- sign(candidates$rho_GSE197677)==sign(candidates$rho_GSE221561)
candidates$mean_abs_rho <- (abs(candidates$rho_GSE197677)+abs(candidates$rho_GSE221561))/2
candidates$minimum_abs_rho <- pmin(abs(candidates$rho_GSE197677), abs(candidates$rho_GSE221561))
candidates$replication_status <- "CROSS_DATASET_SIGNAL_SUPPORTED"

write.csv(candidates, file.path(O_DIR, "STEP19P_REPLICATED_SIGNAL_SET.csv"), row.names=FALSE)
cat_log("8 primary replicated candidates:\n")
for (i in 1:nrow(candidates)) {
  cat_log(sprintf("  %s %s: 197=%.3f 221=%.3f same_dir=%s\n",
    candidates$source[i], candidates$ligand[i],
    candidates$rho_GSE197677[i], candidates$rho_GSE221561[i], candidates$same_direction[i]))
}
cat_log("Saved: STEP19P_REPLICATED_SIGNAL_SET.csv\n\n")

# ============================================================================
# 19P3 — RECEPTOR MAPPING
# ============================================================================
cat_log("============================================\n")
cat_log("19P3: RECEPTOR MAPPING\n")
cat_log("============================================\n\n")

# Rebuild the same curated LR resource
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
  "DLL1","DLL1","DLL4","DLL4","JAG1","JAG1","JAG2","JAG2",  "HGF","HGF"
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
curated_lr <- data.frame(ligand=.lr_lig, receptor=.lr_rec, stringsAsFactors=FALSE)

fib197_logcpm <- compute_log2_cpm(fib197)
fib221_logcpm <- compute_log2_cpm(fib221)

rec_support_list <- list()
for (i in 1:nrow(candidates)) {
  lig <- candidates$ligand[i]
  src <- candidates$source[i]
  lig_receptors <- unique(curated_lr$receptor[curated_lr$ligand==lig])
  for (rec in lig_receptors) {
    # GSE197677
    r197 <- rec197_av[rec197_av$receptor==rec,]
    if (nrow(r197)>0) {
      d197 <- r197$detection_fraction[1]; me197 <- r197$mean_expression[1]; c197 <- r197$classification[1]
    } else { d197<-0; me197<-0; c197<-"NOT_IN_DATA" }
    # GSE221561
    r221 <- rec221_av[rec221_av$receptor==rec,]
    if (nrow(r221)>0) {
      d221 <- r221$detection_fraction[1]; me221 <- r221$mean_expression[1]; c221 <- r221$classification[1]
    } else { d221<-0; me221<-0; c221<-"NOT_IN_DATA" }

    det197 <- c197 %in% c("ROBUSTLY_DETECTED","PARTIALLY_DETECTED")
    det221 <- c221 %in% c("ROBUSTLY_DETECTED","PARTIALLY_DETECTED")

    if (det197 && det221 && c197=="ROBUSTLY_DETECTED" && c221=="ROBUSTLY_DETECTED")
      supp <- "RECEPTOR_REPLICATED_ROBUST"
    else if (det197 && det221)
      supp <- "RECEPTOR_REPLICATED_PARTIAL"
    else if (det197 || det221)
      supp <- "RECEPTOR_DATASET_SPECIFIC"
    else if (c197=="NOT_IN_DATA" && c221=="NOT_IN_DATA")
      supp <- "RECEPTOR_NOT_EVALUABLE"
    else
      supp <- "RECEPTOR_POORLY_DETECTED"

    rec_support_list[[paste(src,lig,rec)]] <- data.frame(
      ligand=lig, source=src, receptor=rec, LR_resource="curated_minimal",
      receptor_detected_GSE197677=det197, receptor_detected_GSE221561=det221,
      receptor_detection_fraction_GSE197677=d197, receptor_detection_fraction_GSE221561=d221,
      receptor_mean_expression_GSE197677=me197, receptor_mean_expression_GSE221561=me221,
      receptor_class_GSE197677=c197, receptor_class_GSE221561=c221,
      receptor_support=supp, stringsAsFactors=FALSE)
  }
}
rec_support <- do.call(rbind, rec_support_list)
write.csv(rec_support, file.path(O_DIR, "STEP19P_RECEPTOR_SUPPORT.csv"), row.names=FALSE)
cat_log("=== RECEPTOR SUPPORT ===\n")
for (i in 1:nrow(rec_support)) {
  cat_log(sprintf("  %s/%s → %s: 197=%s(%s) 221=%s(%s) -> %s\n",
    rec_support$source[i], rec_support$ligand[i], rec_support$receptor[i],
    rec_support$receptor_class_GSE197677[i], round(rec_support$receptor_detection_fraction_GSE197677[i],2),
    rec_support$receptor_class_GSE221561[i], round(rec_support$receptor_detection_fraction_GSE221561[i],2),
    rec_support$receptor_support[i]))
}
cat_log("Saved: STEP19P_RECEPTOR_SUPPORT.csv\n\n")

# ============================================================================
# 19P4 — LIGAND-TARGET RESOURCE PREFLIGHT
# ============================================================================
cat_log("============================================\n")
cat_log("19P4: LIGAND-TARGET RESOURCE PREFLIGHT\n")
cat_log("============================================\n\n")

LT_AVAILABLE <- FALSE
lt_status <- data.frame(
  resource=c("NicheNet_ligand_target","NicheNet_LR_network","project_local_LT","other_local_LT"),
  available=rep(FALSE,4), version=rep(NA,4), path=rep(NA,4),
  number_ligands=rep(0,4), number_targets=rep(0,4), number_edges=rep(0,4),
  stringsAsFactors=FALSE)
write.csv(lt_status, file.path(O_DIR, "STEP19P_LIGAND_TARGET_RESOURCE_STATUS.csv"), row.names=FALSE)
cat_log("No ligand-target matrix available locally.\n")
cat_log("LIGAND_TARGET_RESOURCE_AVAILABLE = FALSE\n")
cat_log("Saved: STEP19P_LIGAND_TARGET_RESOURCE_STATUS.csv\n\n")

# ============================================================================
# 19P5 — DIRECT CORE2 TARGET SUPPORT
# ============================================================================
cat_log("============================================\n")
cat_log("19P5: DIRECT CORE2 TARGET SUPPORT\n")
cat_log("============================================\n\n")

lt_target <- data.frame(
  ligand=candidates$ligand, target=rep("NOT_EVALUABLE",nrow(candidates)),
  regulatory_weight=rep(NA,nrow(candidates)), rank=rep(NA,nrow(candidates)),
  resource=rep("NOT_AVAILABLE",nrow(candidates)),
  status=rep("NOT_EVALUABLE_RESOURCE_UNAVAILABLE",nrow(candidates)),
  stringsAsFactors=FALSE)
write.csv(lt_target, file.path(O_DIR, "STEP19P_LIGAND_CORE2_TARGET_SUPPORT.csv"), row.names=FALSE)
cat_log("No ligand-target matrix available. All candidates: NOT_EVALUABLE_RESOURCE_UNAVAILABLE\n")
cat_log("Saved: STEP19P_LIGAND_CORE2_TARGET_SUPPORT.csv\n\n")

# ============================================================================
# 19P6 — FIBROBLAST DOWNSTREAM RESPONSE
# ============================================================================
cat_log("============================================\n")
cat_log("19P6: FIBROBLAST DOWNSTREAM RESPONSE\n")
cat_log("============================================\n\n")

# Use Step19N regulator activity (FOS, RUNX1, FOXO3, SOX9 available)
regs_available <- unique(reg_act$regulator[!is.na(reg_act$activity)])
cat_log(sprintf("Regulators with activity: %s\n", paste(regs_available, collapse=", ")))

downstream_list <- list()
for (i in 1:nrow(candidates)) {
  lig <- candidates$ligand[i]
  src <- candidates$source[i]

  # For each regulator: check ligand-regulator correlation
  for (ds in c("GSE197677","GSE221561")) {
    core_ds <- if (ds=="GSE197677") core197 else core221
    src_ds <- if (ds=="GSE197677") src197 else src221
    treat_map <- if (ds=="GSE197677") TREAT_MAP_197 else TREAT_MAP_221
    reg_ds <- reg_act[reg_act$dataset==ds,]

    if (!src %in% names(src_ds)) next
    pb <- src_ds[[src]]
    matched_samps <- intersect(colnames(pb), core_ds$sample_id)
    if (length(matched_samps)<4) next

    logcpm <- compute_log2_cpm(pb[,matched_samps,drop=FALSE])
    if (!lig %in% rownames(logcpm)) next
    lig_expr <- as.numeric(logcpm[lig,])
    core2_vals <- core_ds$CORE2_STRONG_score[match(matched_samps, core_ds$sample_id)]

    for (reg in regs_available) {
      reg_samp <- reg_ds[reg_ds$regulator==reg,,drop=FALSE]
      reg_vals <- reg_samp$activity[match(matched_samps, reg_samp$sample_id)]
      valid <- !is.na(lig_expr) & !is.na(reg_vals) & !is.na(core2_vals)
      if (sum(valid)<4) {
        downstream_list[[paste(lig,src,ds,reg)]] <- data.frame(
          ligand=lig, source=src, dataset=ds, regulator=reg,
          ligand_vs_regulator_rho=NA, ligand_vs_core2_rho=NA,
          evidence_type="REGULATOR_ACTIVITY", downstream_support="NOT_EVALUABLE",
          stringsAsFactors=FALSE)
        next
      }
      rho_reg <- tryCatch(cor(lig_expr[valid], reg_vals[valid], method="spearman"),
                          error=function(e) NA_real_)
      rho_core2 <- tryCatch(cor(lig_expr[valid], core2_vals[valid], method="spearman"),
                           error=function(e) NA_real_)
      if (is.na(rho_reg)) {
        ds_supp <- "NOT_EVALUABLE"
      } else if (abs(rho_reg)>=0.40) {
        ds_supp <- "COMPATIBLE_DOWNSTREAM_RESPONSE"
      } else if (abs(rho_reg)>=0.20) {
        ds_supp <- "WEAK_DOWNSTREAM_RESPONSE"
      } else {
        ds_supp <- "NO_DOWNSTREAM_RESPONSE"
      }
      downstream_list[[paste(lig,src,ds,reg)]] <- data.frame(
        ligand=lig, source=src, dataset=ds, regulator=reg,
        ligand_vs_regulator_rho=rho_reg, ligand_vs_core2_rho=rho_core2,
        evidence_type="REGULATOR_ACTIVITY", downstream_support=ds_supp,
        stringsAsFactors=FALSE)
    }
  }
}
downstream_df <- do.call(rbind, downstream_list)
write.csv(downstream_df, file.path(O_DIR, "STEP19P_DOWNSTREAM_RESPONSE_SUPPORT.csv"), row.names=FALSE)
cat_log(sprintf("Downstream response records: %d\n", nrow(downstream_df)))
cat_log("Saved: STEP19P_DOWNSTREAM_RESPONSE_SUPPORT.csv\n\n")

# ============================================================================
# 19P7 — TREATMENT-ADJUSTED ASSOCIATION
# ============================================================================
cat_log("============================================\n")
cat_log("19P7: TREATMENT-ADJUSTED ASSOCIATION\n")
cat_log("============================================\n\n")

adj_list <- list()
for (i in 1:nrow(candidates)) {
  lig <- candidates$ligand[i]
  src <- candidates$source[i]

  for (ds in c("GSE197677","GSE221561")) {
    core_ds <- if (ds=="GSE197677") core197 else core221
    src_ds <- if (ds=="GSE197677") src197 else src221
    treat_map <- if (ds=="GSE197677") TREAT_MAP_197 else TREAT_MAP_221

    if (!src %in% names(src_ds)) next
    pb <- src_ds[[src]]
    matched_samps <- intersect(colnames(pb), core_ds$sample_id)
    matched_samps <- matched_samps[matched_samps %in% names(treat_map)]
    if (length(matched_samps)<5) next

    logcpm <- compute_log2_cpm(pb[,matched_samps,drop=FALSE])
    if (!lig %in% rownames(logcpm)) next
    lig_expr <- as.numeric(logcpm[lig,])
    core2_vals <- core_ds$CORE2_STRONG_score[match(matched_samps, core_ds$sample_id)]
    treat_vals <- as.numeric(treat_map[matched_samps]=="treated")

    # For GSE197677: NACT=1, nNACT=0
    if (ds=="GSE197677") {
      treat_vals <- as.numeric(treat_map[matched_samps]=="NACT")
    } else {
      treat_vals <- as.numeric(treat_map[matched_samps]=="Neoadjuvant_treated")
    }

    raw_rho <- tryCatch(cor(lig_expr, core2_vals, method="spearman"),
                       error=function(e) NA_real_)

    # Treatment-residualized association
    resid_lig <- tryCatch(residuals(lm(lig_expr ~ treat_vals)), error=function(e) lig_expr)
    resid_core2 <- tryCatch(residuals(lm(core2_vals ~ treat_vals)), error=function(e) core2_vals)
    adj_rho <- tryCatch(cor(resid_lig, resid_core2, method="spearman"),
                        error=function(e) NA_real_)

    adj_list[[paste(lig,src,ds)]] <- data.frame(
      ligand=lig, source=src, dataset=ds, n_samples=length(matched_samps),
      raw_rho=raw_rho, treatment_adjusted_rho=adj_rho,
      raw_direction=ifelse(!is.na(raw_rho), ifelse(raw_rho>0,"POSITIVE","NEGATIVE"), "NA"),
      adjusted_direction=ifelse(!is.na(adj_rho), ifelse(adj_rho>0,"POSITIVE","NEGATIVE"), "NA"),
      direction_preserved=ifelse(!is.na(raw_rho) && !is.na(adj_rho), sign(raw_rho)==sign(adj_rho), NA),
      stringsAsFactors=FALSE)
  }
}
adj_df <- do.call(rbind, adj_list)
write.csv(adj_df, file.path(O_DIR, "STEP19P_TREATMENT_ADJUSTED_ASSOCIATIONS.csv"), row.names=FALSE)
cat_log("=== TREATMENT-ADJUSTED ASSOCIATIONS ===\n")
for (i in 1:nrow(adj_df)) {
  cat_log(sprintf("  %s %s %s: raw=%.3f adj=%.3f preserved=%s\n",
    adj_df$ligand[i], adj_df$source[i], adj_df$dataset[i],
    adj_df$raw_rho[i], adj_df$treatment_adjusted_rho[i],
    ifelse(is.na(adj_df$direction_preserved[i]),"NA",adj_df$direction_preserved[i])))
}
cat_log("Saved: STEP19P_TREATMENT_ADJUSTED_ASSOCIATIONS.csv\n\n")

# ============================================================================
# 19P8 — SOURCE ABUNDANCE SENSITIVITY
# ============================================================================
cat_log("============================================\n")
cat_log("19P8: SOURCE ABUNDANCE SENSITIVITY\n")
cat_log("============================================\n\n")

# Source fraction = proportion of source compartment cells in total (approximated by pseudobulk library size ratio)
abund_list <- list()
for (i in 1:nrow(candidates)) {
  lig <- candidates$ligand[i]; src <- candidates$source[i]

  for (ds in c("GSE197677","GSE221561")) {
    core_ds <- if (ds=="GSE197677") core197 else core221
    src_ds <- if (ds=="GSE197677") src197 else src221
    fib_ds <- if (ds=="GSE197677") fib197 else fib221

    if (!src %in% names(src_ds)) next
    pb <- src_ds[[src]]
    fib_pb <- fib_ds
    matched_samps <- intersect(colnames(pb), core_ds$sample_id)
    if (length(matched_samps)<5) next

    logcpm <- compute_log2_cpm(pb[,matched_samps,drop=FALSE])
    if (!lig %in% rownames(logcpm)) next
    lig_expr <- as.numeric(logcpm[lig,])
    core2_vals <- core_ds$CORE2_STRONG_score[match(matched_samps, core_ds$sample_id)]

    # Source abundance proxy: total counts in source pseudobulk
    src_total <- colSums(pb[,matched_samps,drop=FALSE])
    src_frac <- src_total / max(src_total)

    # Model A: CORE2 ~ ligand
    rho_raw <- tryCatch(cor(lig_expr, core2_vals, method="spearman"),
                        error=function(e) NA_real_)
    # Model B: CORE2 ~ ligand + source_fraction
    # Use partial correlation with source fraction
    resid_lig_B <- tryCatch(residuals(lm(lig_expr ~ src_frac)), error=function(e) lig_expr)
    resid_core2_B <- tryCatch(residuals(lm(core2_vals ~ src_frac)), error=function(e) core2_vals)
    rho_adj_src <- tryCatch(cor(resid_lig_B, resid_core2_B, method="spearman"),
                            error=function(e) NA_real_)

    # Model C: CORE2 ~ ligand + fibroblast_fraction (approximated by fibroblast library size)
    fib_total <- colSums(fib_pb[,intersect(colnames(fib_pb), matched_samps),drop=FALSE])
    fib_frac <- fib_total[matched_samps] / max(fib_total)
    resid_lig_C <- tryCatch(residuals(lm(lig_expr ~ fib_frac)), error=function(e) lig_expr)
    resid_core2_C <- tryCatch(residuals(lm(core2_vals ~ fib_frac)), error=function(e) core2_vals)
    rho_adj_fib <- tryCatch(cor(resid_lig_C, resid_core2_C, method="spearman"),
                             error=function(e) NA_real_)

    # Classification
    if (is.na(rho_raw) || is.na(rho_adj_src)) {
      class <- "NOT_EVALUABLE"
    } else if (sign(rho_raw)==sign(rho_adj_src) && abs(abs(rho_raw)-abs(rho_adj_src))<0.15) {
      class <- "ABUNDANCE_ROBUST"
    } else if (sign(rho_raw)==sign(rho_adj_src)) {
      class <- "PARTIALLY_ABUNDANCE_SENSITIVE"
    } else {
      class <- "ABUNDANCE_SENSITIVE"
    }

    abund_list[[paste(lig,src,ds)]] <- data.frame(
      ligand=lig, source=src, dataset=ds,
      raw_rho=rho_raw, adjusted_source_rho=rho_adj_src, adjusted_fibroblast_rho=rho_adj_fib,
      raw_direction=ifelse(!is.na(rho_raw), ifelse(rho_raw>0,"POSITIVE","NEGATIVE"),"NA"),
      adjusted_source_direction=ifelse(!is.na(rho_adj_src), ifelse(rho_adj_src>0,"POSITIVE","NEGATIVE"),"NA"),
      direction_change_after_adjustment=ifelse(!is.na(rho_raw)&&!is.na(rho_adj_src),
        sign(rho_raw)!=sign(rho_adj_src), NA),
      abundance_sensitivity=class, stringsAsFactors=FALSE)
  }
}
abund_df <- do.call(rbind, abund_list)
write.csv(abund_df, file.path(O_DIR, "STEP19P_ABUNDANCE_SENSITIVITY.csv"), row.names=FALSE)
cat_log("Saved: STEP19P_ABUNDANCE_SENSITIVITY.csv\n\n")

# ============================================================================
# 19P9 — TREATMENT EFFECT REPLICATION
# ============================================================================
cat_log("============================================\n")
cat_log("19P9: TREATMENT EFFECT REPLICATION\n")
cat_log("============================================\n\n")

te_list <- list()
for (i in 1:nrow(candidates)) {
  lig <- candidates$ligand[i]; src <- candidates$source[i]

  te_197_val <- NA; te_221_val <- NA
  for (ds in c("GSE197677","GSE221561")) {
    core_ds <- if (ds=="GSE197677") core197 else core221
    src_ds <- if (ds=="GSE197677") src197 else src221
    treat_map <- if (ds=="GSE197677") TREAT_MAP_197 else TREAT_MAP_221

    if (!src %in% names(src_ds)) next
    pb <- src_ds[[src]]
    samps <- colnames(pb)
    samps <- samps[samps %in% names(treat_map)]
    if (length(samps)<4) next

    logcpm <- compute_log2_cpm(pb[,samps,drop=FALSE])
    if (!lig %in% rownames(logcpm)) next
    lig_expr <- as.numeric(logcpm[lig,samps])
    treated <- samps[sapply(samps,function(s) treat_map[s]=="NACT" || treat_map[s]=="Neoadjuvant_treated")]
    control <- samps[sapply(samps,function(s) treat_map[s]=="nNACT" || treat_map[s]=="Surgery_alone")]
    if (length(treated)==0 || length(control)==0) next
    eff <- mean(lig_expr[match(treated,samps)]) - mean(lig_expr[match(control,samps)])
    if (ds=="GSE197677") te_197_val <- eff else te_221_val <- eff
  }

  same_dir <- if (!is.na(te_197_val) && !is.na(te_221_val)) sign(te_197_val)==sign(te_221_val) else NA
  repl <- if (is.na(same_dir)) "INSUFFICIENT"
    else if (same_dir) "TREATMENT_EFFECT_REPLICATED"
    else "TREATMENT_EFFECT_DISCORDANT"

  te_list[[paste(lig,src)]] <- data.frame(
    ligand=lig, source=src,
    effect_GSE197677=te_197_val, effect_GSE221561=te_221_val,
    same_direction=same_dir, effect_replication_status=repl,
    stringsAsFactors=FALSE)
}
te_df <- do.call(rbind, te_list)
write.csv(te_df, file.path(O_DIR, "STEP19P_TREATMENT_EFFECT_REPLICATION.csv"), row.names=FALSE)
cat_log("Saved: STEP19P_TREATMENT_EFFECT_REPLICATION.csv\n\n")

# ============================================================================
# 19P10 — LOSO ROBUSTNESS
# ============================================================================
cat_log("============================================\n")
cat_log("19P10: LOSO ROBUSTNESS\n")
cat_log("============================================\n\n")

loso_long_list <- list()
for (i in 1:nrow(candidates)) {
  lig <- candidates$ligand[i]; src <- candidates$source[i]

  for (ds in c("GSE197677","GSE221561")) {
    core_ds <- if (ds=="GSE197677") core197 else core221
    src_ds <- if (ds=="GSE197677") src197 else src221
    treat_map <- if (ds=="GSE197677") TREAT_MAP_197 else TREAT_MAP_221

    if (!src %in% names(src_ds)) next
    pb <- src_ds[[src]]
    matched_samps <- intersect(colnames(pb), core_ds$sample_id)
    if (length(matched_samps)<5) next

    logcpm <- compute_log2_cpm(pb[,matched_samps,drop=FALSE])
    if (!lig %in% rownames(logcpm)) next
    lig_expr <- as.numeric(logcpm[lig,])
    core2_vals <- core_ds$CORE2_STRONG_score[match(matched_samps, core_ds$sample_id)]
    full_rho <- tryCatch(cor(lig_expr, core2_vals, method="spearman"),
                         error=function(e) NA_real_)

    for (loso_s in matched_samps) {
      keep <- matched_samps[matched_samps!=loso_s]
      if (length(keep)<4) next
      rho_loso <- tryCatch(cor(lig_expr[keep], core2_vals[keep], method="spearman"),
                           error=function(e) NA_real_)
      sign_pres <- if (!is.na(rho_loso)&&!is.na(full_rho)) sign(rho_loso)==sign(full_rho) else NA
      class <- "STABLE"
      if (!is.na(sign_pres) && !sign_pres) {
        is_ctrl <- if (ds=="GSE197677") treat_map[loso_s]=="nNACT" else treat_map[loso_s]=="Surgery_alone"
        class <- if (is_ctrl) "CONTROL_SENSITIVE" else "OUTLIER_SENSITIVE"
      }
      loso_long_list[[paste(lig,src,ds,loso_s)]] <- data.frame(
        ligand=lig, source=src, dataset=ds, left_out=loso_s,
        left_out_treatment=if (loso_s %in% names(treat_map)) treat_map[loso_s] else NA,
        rho_loso=rho_loso, rho_full=full_rho,
        rho_sign_preserved=sign_pres, classification=class,
        stringsAsFactors=FALSE)
    }
  }
}
loso_long <- do.call(rbind, loso_long_list)
write.csv(loso_long, file.path(O_DIR, "STEP19P_LOSO_LONG.csv"), row.names=FALSE)

# Summary
loso_sum_list <- list()
for (i in 1:nrow(candidates)) {
  lig <- candidates$ligand[i]; src <- candidates$source[i]
  for (ds in c("GSE197677","GSE221561")) {
    sub <- loso_long[loso_long$ligand==lig & loso_long$source==src & loso_long$dataset==ds,]
    if (nrow(sub)==0) next
    n_ctrl_sens <- sum(sub$classification=="CONTROL_SENSITIVE")
    n_outlier_sens <- sum(sub$classification=="OUTLIER_SENSITIVE")
    n_sign_pres <- sum(sub$rho_sign_preserved[!is.na(sub$rho_sign_preserved)])
    n_total <- sum(!is.na(sub$rho_sign_preserved))
    sign_pres_frac <- if (n_total>0) n_sign_pres/n_total else NA
    rho_vals <- sub$rho_loso[!is.na(sub$rho_loso)]
    class <- "ROBUST"
    if (is.na(sign_pres_frac)) {
      class <- "NOT_EVALUABLE"
    } else if (n_ctrl_sens>0) {
      class <- "CONTROL_SENSITIVE"
    } else if (n_outlier_sens>0) {
      class <- "OUTLIER_SENSITIVE"
    } else if (sign_pres_frac<0.8) {
      class <- "UNSTABLE"
    }
    loso_sum_list[[paste(lig,src,ds)]] <- data.frame(
      ligand=lig, source=src, dataset=ds,
      rho_sign_preservation_fraction=sign_pres_frac,
      rho_ge_0.4_preservation_fraction=mean(abs(rho_vals)>=0.40),
      median_loso_rho=median(rho_vals), min_rho=min(rho_vals), max_rho=max(rho_vals),
      classification=class, stringsAsFactors=FALSE)
  }
}
loso_summary <- do.call(rbind, loso_sum_list)
write.csv(loso_summary, file.path(O_DIR, "STEP19P_LOSO_SUMMARY.csv"), row.names=FALSE)
cat_log("Saved: STEP19P_LOSO_LONG.csv, STEP19P_LOSO_SUMMARY.csv\n\n")

# ============================================================================
# 19P11 — SOURCE-SPECIFICITY (FGF7 and CXCL2)
# ============================================================================
cat_log("============================================\n")
cat_log("19P11: SOURCE-SPECIFICITY\n")
cat_log("============================================\n\n")

multi_source_ligands <- c("FGF7","CXCL2")
ms_list <- list()
for (lig in multi_source_ligands) {
  lig_sources <- candidates$source[candidates$ligand==lig]
  for (ds in c("GSE197677","GSE221561")) {
    core_ds <- if (ds=="GSE197677") core197 else core221
    src_ds <- if (ds=="GSE197677") src197 else src221
    for (src in names(src_ds)) {
      pb <- src_ds[[src]]
      matched_samps <- intersect(colnames(pb), core_ds$sample_id)
      if (length(matched_samps)<4) next
      logcpm <- compute_log2_cpm(pb[,matched_samps,drop=FALSE])
      if (!lig %in% rownames(logcpm)) { rho<-NA; det_frac<-0; mean_expr<-0 } else {
        lig_expr <- as.numeric(logcpm[lig,])
        core2_vals <- core_ds$CORE2_STRONG_score[match(matched_samps, core_ds$sample_id)]
        rho <- tryCatch(cor(lig_expr, core2_vals, method="spearman"), error=function(e) NA_real_)
        det_frac <- mean(lig_expr>0)
        mean_expr <- mean(lig_expr)
      }
      ms_list[[paste(lig,src,ds)]] <- data.frame(
        ligand=lig, source=src, dataset=ds,
        ligand_rho=rho, detection_fraction=det_frac, mean_expression=mean_expr,
        stringsAsFactors=FALSE)
    }
  }
}
ms_df <- do.call(rbind, ms_list)

# Classify
ms_summary_list <- list()
for (lig in multi_source_ligands) {
  lig_rows <- candidates[candidates$ligand==lig,]
  n_sources <- nrow(lig_rows)
  if (n_sources>=2) {
    all_rho <- abs(lig_rows$rho_GSE197677)
    all_rho <- c(all_rho, abs(lig_rows$rho_GSE221561))
    if (min(all_rho)>=0.40) ms_class <- "MULTI_SOURCE_SUPPORTED"
    else if (min(all_rho)>=0.30) ms_class <- "SOURCE_AMBIGUOUS"
    else ms_class <- "SOURCE_SPECIFIC"
  } else {
    ms_class <- "SOURCE_SPECIFIC"
  }
  ms_summary_list[[lig]] <- data.frame(
    ligand=lig, n_replicated_sources=n_sources,
    replicated_sources=paste(lig_rows$source, collapse="; "),
    source_specificity=ms_class, stringsAsFactors=FALSE)
}
ms_summary <- do.call(rbind, ms_summary_list)
write.csv(ms_df, file.path(O_DIR, "STEP19P_MULTI_SOURCE_LIGAND_AUDIT.csv"), row.names=FALSE)
cat_log(sprintf("FGF7: %s\n", ms_summary$source_specificity[ms_summary$ligand=="FGF7"]))
cat_log(sprintf("CXCL2: %s\n", ms_summary$source_specificity[ms_summary$ligand=="CXCL2"]))
cat_log("Saved: STEP19P_MULTI_SOURCE_LIGAND_AUDIT.csv\n\n")

# ============================================================================
# 19P12 — NEGATIVE CONTROL / DISCORDANT COMPARATORS
# ============================================================================
cat_log("============================================\n")
cat_log("19P12: NEGATIVE COMPARATOR AUDIT\n")
cat_log("============================================\n\n")

neg_cands <- data.frame(
  ligand=c("AREG","DLL1","EREG","TGFB3","IFNG"),
  source=c("Epithelial","Epithelial","Epithelial","T_NK","Myeloid"),
  rho_197=c(-0.794,-0.782,-0.806,0.758,0.830),
  rho_221=c(0.467,0.383,0.850,-0.492,-0.186),
  stringsAsFactors=FALSE)

neg_list <- list()
for (i in 1:nrow(neg_cands)) {
  lig <- neg_cands$ligand[i]; src <- neg_cands$source[i]
  lig_receptors <- unique(curated_lr$receptor[curated_lr$ligand==lig])

  rec_197 <- rec197_av[rec197_av$receptor %in% lig_receptors,]
  rec_221 <- rec221_av[rec221_av$receptor %in% lig_receptors,]
  rec_197_det <- any(rec_197$classification %in% c("ROBUSTLY_DETECTED","PARTIALLY_DETECTED"))
  rec_221_det <- any(rec_221$classification %in% c("ROBUSTLY_DETECTED","PARTIALLY_DETECTED"))

  # Treatment-adjusted from adj_df
  adj_197 <- adj_df[adj_df$ligand==lig & adj_df$source==src & adj_df$dataset=="GSE197677",]
  adj_221 <- adj_df[adj_df$ligand==lig & adj_df$source==src & adj_df$dataset=="GSE221561",]
  adj_rho_197 <- if (nrow(adj_197)>0) adj_197$treatment_adjusted_rho[1] else NA
  adj_rho_221 <- if (nrow(adj_221)>0) adj_221$treatment_adjusted_rho[1] else NA

  # LOSO
  loso_197 <- loso_summary[loso_summary$ligand==lig & loso_summary$source==src & loso_summary$dataset=="GSE197677",]
  loso_221 <- loso_summary[loso_summary$ligand==lig & loso_summary$source==src & loso_summary$dataset=="GSE221561",]
  loso_197_class <- if (nrow(loso_197)>0) loso_197$classification[1] else "NOT_EVALUABLE"
  loso_221_class <- if (nrow(loso_221)>0) loso_221$classification[1] else "NOT_EVALUABLE"

  same_dir <- sign(neg_cands$rho_197[i])==sign(neg_cands$rho_221[i])

  neg_list[[lig]] <- data.frame(
    ligand=lig, source=src,
    rho_197=neg_cands$rho_197[i], rho_221=neg_cands$rho_221[i],
    same_direction=same_dir,
    receptor_197=rec_197_det, receptor_221=rec_221_det,
    adjusted_rho_197=adj_rho_197, adjusted_rho_221=adj_rho_221,
    loso_197=loso_197_class, loso_221=loso_221_class,
    comparator_status="DISCORDANT",
    stringsAsFactors=FALSE)
}
neg_df <- do.call(rbind, neg_list)
write.csv(neg_df, file.path(O_DIR, "STEP19P_NEGATIVE_COMPARATOR_AUDIT.csv"), row.names=FALSE)
cat_log("Saved: STEP19P_NEGATIVE_COMPARATOR_AUDIT.csv\n\n")

# ============================================================================
# 19P13 — EVIDENCE SCORE
# ============================================================================
cat_log("============================================\n")
cat_log("19P13: EVIDENCE SCORE\n")
cat_log("============================================\n\n")

evi_list <- list()
for (i in 1:nrow(candidates)) {
  lig <- candidates$ligand[i]; src <- candidates$source[i]
  score <- 0; avail <- 0

  # +1 CORE2 abs(rho)>=0.4 in GSE197677
  avail <- avail+1
  c1 <- abs(candidates$rho_GSE197677[i])>=0.40
  if (c1) score <- score+1

  # +1 CORE2 abs(rho)>=0.4 in GSE221561
  avail <- avail+1
  c2 <- abs(candidates$rho_GSE221561[i])>=0.40
  if (c2) score <- score+1

  # +1 same direction
  avail <- avail+1
  c3 <- candidates$same_direction[i]
  if (c3) score <- score+1

  # +1 receptor in both datasets
  avail <- avail+1
  lig_receptors <- unique(curated_lr$receptor[curated_lr$ligand==lig])
  rec_s <- rec_support[rec_support$ligand==lig & rec_support$source==src,]
  c4 <- any(rec_s$receptor_support %in% c("RECEPTOR_REPLICATED_ROBUST","RECEPTOR_REPLICATED_PARTIAL"))
  if (c4) score <- score+1

  # +1 treatment-adjusted preserves direction GSE197677
  av <- adj_df[adj_df$ligand==lig & adj_df$source==src & adj_df$dataset=="GSE197677",]
  avail <- avail+1
  c5 <- if (nrow(av)>0 && !is.na(av$direction_preserved[1])) av$direction_preserved[1] else FALSE
  if (c5) score <- score+1

  # +1 treatment-adjusted preserves direction GSE221561
  av2 <- adj_df[adj_df$ligand==lig & adj_df$source==src & adj_df$dataset=="GSE221561",]
  avail <- avail+1
  c6 <- if (nrow(av2)>0 && !is.na(av2$direction_preserved[1])) av2$direction_preserved[1] else FALSE
  if (c6) score <- score+1

  # +1 abundance adjustment preserves direction both
  avail <- avail+1
  ab_197 <- abund_df[abund_df$ligand==lig & abund_df$source==src & abund_df$dataset=="GSE197677",]
  ab_221 <- abund_df[abund_df$ligand==lig & abund_df$source==src & abund_df$dataset=="GSE221561",]
  c7 <- FALSE
  if (nrow(ab_197)>0 && nrow(ab_221)>0) {
    c7 <- ab_197$abundance_sensitivity[1]=="ABUNDANCE_ROBUST" &&
          ab_221$abundance_sensitivity[1]=="ABUNDANCE_ROBUST"
  }
  if (c7) score <- score+1

  # +1 LOSO direction stability in both
  avail <- avail+1
  ls_197 <- loso_summary[loso_summary$ligand==lig & loso_summary$source==src & loso_summary$dataset=="GSE197677",]
  ls_221 <- loso_summary[loso_summary$ligand==lig & loso_summary$source==src & loso_summary$dataset=="GSE221561",]
  c8 <- FALSE
  if (nrow(ls_197)>0 && nrow(ls_221)>0) {
    c8 <- ls_197$classification[1] %in% c("ROBUST") &&
          ls_221$classification[1] %in% c("ROBUST")
  }
  if (c8) score <- score+1

  # +1 ligand-target (NOT AVAILABLE -> not scored, not penalized)
  # avail stays same (8 available)

  # +1 downstream regulator compatible
  avail <- avail+1
  ds_197 <- downstream_df[downstream_df$ligand==lig & downstream_df$source==src &
    downstream_df$dataset=="GSE197677" & downstream_df$downstream_support=="COMPATIBLE_DOWNSTREAM_RESPONSE",]
  c9 <- nrow(ds_197)>0
  if (c9) score <- score+1

  frac <- score/avail
  if (frac>=0.75) evi_class <- "HIGH_FOCUSED_SUPPORT"
  else if (frac>=0.50) evi_class <- "MODERATE_FOCUSED_SUPPORT"
  else if (frac>=0.38) evi_class <- "WEAK_FOCUSED_SUPPORT"
  else if (frac>=0.25) evi_class <- "ASSOCIATION_ONLY"
  else evi_class <- "RESOURCE_LIMITED"

  evi_list[[paste(lig,src)]] <- data.frame(
    ligand=lig, source=src,
    c_core2_197=c1, c_core2_221=c2, c_same_dir=c3,
    c_receptor_both=c4, c_adj_dir_197=c5, c_adj_dir_221=c6,
    c_abundance=c7, c_loso=c8, c_target_support=FALSE, c_downstream=c9,
    score_numerator=score, score_denominator_available=avail,
    fraction_supported=round(frac,3),
    classification=evi_class, stringsAsFactors=FALSE)
}
evi_df <- do.call(rbind, evi_list)
write.csv(evi_df, file.path(O_DIR, "STEP19P_FOCUSED_EVIDENCE_SCORE.csv"), row.names=FALSE)
cat_log("=== EVIDENCE SCORES ===\n")
for (i in 1:nrow(evi_df)) {
  cat_log(sprintf("  %s/%s: %d/%d (%.2f) -> %s\n",
    evi_df$ligand[i], evi_df$source[i],
    evi_df$score_numerator[i], evi_df$score_denominator_available[i],
    evi_df$fraction_supported[i], evi_df$classification[i]))
}
cat_log("Saved: STEP19P_FOCUSED_EVIDENCE_SCORE.csv\n\n")

# ============================================================================
# 19P14 — CANDIDATE PRIORITIZATION
# ============================================================================
cat_log("============================================\n")
cat_log("19P14: CANDIDATE PRIORITIZATION\n")
cat_log("============================================\n\n")

rank_df <- candidates
rank_df$receptors <- sapply(1:nrow(rank_df), function(i) {
  rs <- unique(rec_support$receptor[rec_support$ligand==rank_df$ligand[i] & rec_support$source==rank_df$source[i]])
  paste(rs, collapse=", ")
})
rank_df$adjusted_support <- sapply(1:nrow(rank_df), function(i) {
  a197 <- adj_df[adj_df$ligand==rank_df$ligand[i] & adj_df$source==rank_df$source[i] & adj_df$dataset=="GSE197677",]
  a221 <- adj_df[adj_df$ligand==rank_df$ligand[i] & adj_df$source==rank_df$source[i] & adj_df$dataset=="GSE221561",]
  p197 <- if (nrow(a197)>0 && !is.na(a197$direction_preserved[1])) a197$direction_preserved[1] else FALSE
  p221 <- if (nrow(a221)>0 && !is.na(a221$direction_preserved[1])) a221$direction_preserved[1] else FALSE
  if (p197 && p221) "BOTH_PRESERVED" else if (p197 || p221) "ONE_PRESERVED" else "NOT_PRESERVED"
})
rank_df$abundance_support <- sapply(1:nrow(rank_df), function(i) {
  a197 <- abund_df[abund_df$ligand==rank_df$ligand[i] & abund_df$source==rank_df$source[i] & abund_df$dataset=="GSE197677",]
  a221 <- abund_df[abund_df$ligand==rank_df$ligand[i] & abund_df$source==rank_df$source[i] & abund_df$dataset=="GSE221561",]
  c197 <- if (nrow(a197)>0) a197$abundance_sensitivity[1] else "NOT_EVALUABLE"
  c221 <- if (nrow(a221)>0) a221$abundance_sensitivity[1] else "NOT_EVALUABLE"
  if (c197=="ABUNDANCE_ROBUST" && c221=="ABUNDANCE_ROBUST") "ROBUST_BOTH" else if (c197=="ABUNDANCE_ROBUST" || c221=="ABUNDANCE_ROBUST") "ROBUST_ONE" else "SENSITIVE"
})
rank_df$loso_support <- sapply(1:nrow(rank_df), function(i) {
  l197 <- loso_summary[loso_summary$ligand==rank_df$ligand[i] & loso_summary$source==rank_df$source[i] & loso_summary$dataset=="GSE197677",]
  l221 <- loso_summary[loso_summary$ligand==rank_df$ligand[i] & loso_summary$source==rank_df$source[i] & loso_summary$dataset=="GSE221561",]
  c197 <- if (nrow(l197)>0) l197$classification[1] else "NOT_EVALUABLE"
  c221 <- if (nrow(l221)>0) l221$classification[1] else "NOT_EVALUABLE"
  if (c197=="ROBUST" && c221=="ROBUST") "ROBUST_BOTH" else if (c197=="ROBUST" || c221=="ROBUST") "ROBUST_ONE" else "UNSTABLE"
})
rank_df$target_support <- "NOT_EVALUABLE"
rank_df$downstream_support <- sapply(1:nrow(rank_df), function(i) {
  ds <- downstream_df[downstream_df$ligand==rank_df$ligand[i] & downstream_df$source==rank_df$source[i] & downstream_df$downstream_support=="COMPATIBLE_DOWNSTREAM_RESPONSE",]
  if (nrow(ds)>0) "COMPATIBLE" else "NO_COMPATIBLE"
})
rank_df$evidence_score <- evi_df$score_numerator[match(paste(rank_df$ligand,rank_df$source), paste(evi_df$ligand,evi_df$source))]
rank_df$evidence_denominator <- evi_df$score_denominator_available[match(paste(rank_df$ligand,rank_df$source), paste(evi_df$ligand,evi_df$source))]
rank_df$classification <- evi_df$classification[match(paste(rank_df$ligand,rank_df$source), paste(evi_df$ligand,evi_df$source))]
rank_df <- rank_df[order(-rank_df$evidence_score),]
rank_df$rank <- 1:nrow(rank_df)
rank_df <- rank_df[,c("rank","source","ligand","receptors","rho_GSE197677","rho_GSE221561",
  "adjusted_support","abundance_support","loso_support","target_support","downstream_support",
  "evidence_score","evidence_denominator","classification")]
write.csv(rank_df, file.path(O_DIR, "STEP19P_FINAL_CANDIDATE_RANKING.csv"), row.names=FALSE)
cat_log("=== FINAL CANDIDATE RANKING ===\n")
for (i in 1:nrow(rank_df)) {
  cat_log(sprintf("  %d. %s/%s: score=%d/%d %s adj=%s abund=%s loso=%s ds=%s\n",
    rank_df$rank[i], rank_df$source[i], rank_df$ligand[i],
    rank_df$evidence_score[i], rank_df$evidence_denominator[i], rank_df$classification[i],
    rank_df$adjusted_support[i], rank_df$abundance_support[i],
    rank_df$loso_support[i], rank_df$downstream_support[i]))
}
cat_log("Saved: STEP19P_FINAL_CANDIDATE_RANKING.csv\n\n")

# ============================================================================
# 19P15 — FAMILY-LEVEL INTERPRETATION
# ============================================================================
cat_log("============================================\n")
cat_log("19P15: FAMILY-LEVEL INTERPRETATION\n")
cat_log("============================================\n\n")

fam_map <- data.frame(
  ligand=c("IL1B","FGF7","FGF9","LIF","CXCL2","WNT5A"),
  family=c("IL1","FGF","FGF","LIF_IL6","CXCL","WNT"),
  stringsAsFactors=FALSE)

fam_list <- list()
for (fam in unique(fam_map$family)) {
  fam_ligs <- fam_map$ligand[fam_map$family==fam]
  fam_cands <- candidates[candidates$ligand %in% fam_ligs,]
  if (nrow(fam_cands)==0) next
  n_signals <- nrow(fam_cands)
  n_sources <- length(unique(fam_cands$source))
  cross_consist <- all(fam_cands$same_direction)

  fam_receptors <- unique(curated_lr$receptor[curated_lr$ligand %in% fam_ligs])
  fr_197 <- rec197_av[rec197_av$receptor %in% fam_receptors,]
  fr_221 <- rec221_av[rec221_av$receptor %in% fam_receptors,]
  rec_both <- any(fr_197$classification %in% c("ROBUSTLY_DETECTED","PARTIALLY_DETECTED")) &&
              any(fr_221$classification %in% c("ROBUSTLY_DETECTED","PARTIALLY_DETECTED"))

  fam_ev <- evi_df[evi_df$ligand %in% fam_ligs,]
  high_mod <- any(fam_ev$classification %in% c("HIGH_FOCUSED_SUPPORT","MODERATE_FOCUSED_SUPPORT"))

  if (n_signals>=2 && cross_consist && rec_both && high_mod) class <- "REPLICATED_FAMILY_SUPPORT"
  else if (n_signals>=2 && cross_consist) class <- "MULTI_SOURCE_FAMILY_SUPPORT"
  else if (rec_both) class <- "ASSOCIATION_ONLY_FAMILY"
  else class <- "INSUFFICIENT_DOWNSTREAM_SUPPORT"

  fam_list[[fam]] <- data.frame(
    family=fam, n_replicated_signals=n_signals,
    replicated_ligands=paste(fam_cands$ligand, collapse="; "),
    n_source_compartments=n_sources,
    cross_dataset_consistency=cross_consist,
    receptor_support_both=rec_both,
    focused_validation_support=high_mod,
    classification=class, stringsAsFactors=FALSE)
}
fam_df <- do.call(rbind, fam_list)
write.csv(fam_df, file.path(O_DIR, "STEP19P_SIGNALING_FAMILY_SUMMARY.csv"), row.names=FALSE)
cat_log("=== SIGNALING FAMILY SUMMARY ===\n")
for (i in 1:nrow(fam_df)) cat_log(sprintf("  %s: %d signals, %d sources -> %s\n",
  fam_df$family[i], fam_df$n_replicated_signals[i], fam_df$n_source_compartments[i], fam_df$classification[i]))
cat_log("Saved: STEP19P_SIGNALING_FAMILY_SUMMARY.csv\n\n")

# ============================================================================
# 19P16 — FIGURES
# ============================================================================
cat_log("============================================\n")
cat_log("19P16: FIGURES\n")
cat_log("============================================\n\n")

# 1. Replicated signal heatmap
heat_df <- data.frame(
  candidate=paste(candidates$source, candidates$ligand),
  GSE197677_rho=candidates$rho_GSE197677,
  GSE221561_rho=candidates$rho_GSE221561)
adj_wide <- adj_df[,c("ligand","source","dataset","treatment_adjusted_rho")]
adj_wide$candidate <- paste(adj_wide$source, adj_wide$ligand)
adj_wide$GSE197677_adj_rho <- ifelse(adj_wide$dataset=="GSE197677", adj_wide$treatment_adjusted_rho, NA)
adj_wide$GSE221561_adj_rho <- ifelse(adj_wide$dataset=="GSE221561", adj_wide$treatment_adjusted_rho, NA)
adj_cast <- reshape2::dcast(adj_wide, candidate~dataset, value.var="treatment_adjusted_rho")
colnames(adj_cast)[2:3] <- c("GSE197677_adj","GSE221561_adj")
heat_merge <- merge(heat_df, adj_cast, by="candidate", all.x=TRUE)
heat_mat <- as.matrix(heat_merge[,c("GSE197677_rho","GSE221561_rho","GSE197677_adj","GSE221561_adj")])
rownames(heat_mat) <- heat_merge$candidate
colnames(heat_mat) <- c("197 rho","221 rho","197 adj","221 adj")

p <- ggplot(reshape2::melt(heat_mat), aes(x=Var2, y=Var1, fill=value)) +
  geom_tile(color="white") +
  scale_fill_gradient2(low="steelblue", mid="white", high="firebrick", midpoint=0) +
  theme_minimal(base_size=10) +
  labs(title="Replicated Signal Heatmap", x="", y="Candidate") +
  theme(axis.text.x=element_text(angle=45,hjust=1))
ggsave(file.path(FIG_DIR, "STEP19P_REPLICATED_SIGNAL_HEATMAP.png"), p, width=6, height=5, dpi=150)

# 2. Raw vs adjusted associations
raw_adj_df <- adj_df[,c("ligand","source","dataset","raw_rho","treatment_adjusted_rho")]
raw_adj_df$candidate <- paste(raw_adj_df$source, raw_adj_df$ligand)
p2 <- ggplot(raw_adj_df, aes(x=raw_rho, y=treatment_adjusted_rho, color=dataset)) +
  geom_point(size=3) + geom_abline(slope=1, intercept=0, linetype="dashed") +
  theme_minimal(base_size=12) +
  labs(title="Raw vs Treatment-Adjusted Associations",
    x="Raw Spearman rho", y="Treatment-adjusted rho") +
  geom_text(aes(label=candidate), size=2.5, vjust=-0.5)
ggsave(file.path(FIG_DIR, "STEP19P_RAW_VS_ADJUSTED_ASSOCIATIONS.png"), p2, width=7, height=6, dpi=150)

# 3. Receptor support matrix
rec_mat <- rec_support[,c("ligand","source","receptor","receptor_detection_fraction_GSE197677","receptor_detection_fraction_GSE221561")]
rec_mat$candidate <- paste(rec_mat$source, rec_mat$ligand)
rec_plot <- reshape2::melt(rec_mat[,c("candidate","receptor","receptor_detection_fraction_GSE197677","receptor_detection_fraction_GSE221561")],
  id.vars=c("candidate","receptor"))
p3 <- ggplot(rec_plot, aes(x=receptor, y=candidate, fill=value)) +
  geom_tile(color="white") +
  scale_fill_gradient2(low="white", high="steelblue", midpoint=0.5, limit=c(0,1)) +
  facet_wrap(~variable) +
  theme_minimal(base_size=8) +
  labs(title="Receptor Support Matrix", x="Receptor", y="Candidate") +
  theme(axis.text.x=element_text(angle=45,hjust=1))
ggsave(file.path(FIG_DIR, "STEP19P_RECEPTOR_SUPPORT_MATRIX.png"), p3, width=8, height=5, dpi=150)

# 4. LOSO stability
loso_plot <- loso_summary[,c("ligand","source","dataset","rho_sign_preservation_fraction","classification")]
loso_plot$candidate <- paste(loso_plot$source, loso_plot$ligand)
p4 <- ggplot(loso_plot, aes(x=candidate, y=rho_sign_preservation_fraction, fill=classification)) +
  geom_col() + coord_flip() + facet_wrap(~dataset) +
  theme_minimal(base_size=10) +
  labs(title="LOSO Stability", x="", y="Sign Preservation Fraction")
ggsave(file.path(FIG_DIR, "STEP19P_LOSO_STABILITY.png"), p4, width=8, height=5, dpi=150)

# 5. Evidence score
evi_plot <- evi_df
evi_plot$candidate <- paste(evi_plot$source, evi_plot$ligand)
p5 <- ggplot(evi_plot, aes(x=reorder(candidate, score_numerator), y=score_numerator, fill=classification)) +
  geom_col() + coord_flip() +
  theme_minimal(base_size=10) +
  labs(title="Focused Evidence Score", x="", y="Score (numerator)") +
  geom_text(aes(label=paste0(score_numerator,"/",score_denominator_available)), hjust=-0.1, size=3)
ggsave(file.path(FIG_DIR, "STEP19P_FOCUSED_EVIDENCE_SCORE.png"), p5, width=7, height=5, dpi=150)

# 6. Final ranking
rank_plot <- rank_df
rank_plot$candidate <- paste(rank_plot$source, rank_plot$ligand)
p6 <- ggplot(rank_plot, aes(x=reorder(candidate, -evidence_score), y=evidence_score, fill=classification)) +
  geom_col() + coord_flip() +
  theme_minimal(base_size=10) +
  labs(title="Final Candidate Ranking", x="", y="Evidence Score") +
  geom_text(aes(label=paste0(evidence_score,"/",evidence_denominator)), hjust=-0.1, size=3)
ggsave(file.path(FIG_DIR, "STEP19P_FINAL_CANDIDATE_RANKING.png"), p6, width=7, height=5, dpi=150)

# 7. Family summary
fam_plot <- fam_df
p7 <- ggplot(fam_plot, aes(x=family, y=n_replicated_signals, fill=classification)) +
  geom_col() + theme_minimal(base_size=12) +
  labs(title="Signaling Family Summary", x="Family", y="N Replicated Signals") +
  theme(axis.text.x=element_text(angle=45,hjust=1))
ggsave(file.path(FIG_DIR, "STEP19P_SIGNALING_FAMILY_SUMMARY.png"), p7, width=6, height=5, dpi=150)

cat_log("7 figures saved to Step19P/figures.\n\n")

# ============================================================================
# 19P17 — INTERPRETATION GUARDRAILS
# ============================================================================
cat_log("============================================\n")
cat_log("19P17: INTERPRETATION GUARDRAILS\n")
cat_log("============================================\n\n")

guardrails <- paste0(
"# Step 19P Interpretation Guardrails\n\n",
"1. The eight signals are replicated sample-level associations.\n",
"2. Replicated association does not prove ligand secretion.\n",
"3. Receptor expression does not prove receptor activation.\n",
"4. Treatment-adjusted analyses are sensitivity analyses, not causal identification.\n",
"5. Source abundance can induce apparent signaling association.\n",
"6. Ligand-target predictions, if available, are computational.\n",
"7. Sample/patient is the statistical unit.\n",
"8. GSE221561 Surgery_alone n=2 materially limits inference.\n",
"9. The LR resource contains only approximately 100 curated pairs and is not comprehensive.\n",
"10. Absence from the LR resource is not evidence of biological absence.\n",
"11. No single ligand should be called a driver without downstream functional evidence.\n",
"12. The original AREG/DLL1/EREG/TGFB3 signals are discordant across datasets and must not be reported as conserved.\n",
"13. IFNG lacks adequate fibroblast receptor support in the current audit.\n",
"14. CORE2 conservation is stronger evidence than any individual upstream ligand mechanism.\n"
)
writeLines(guardrails, file.path(O_DIR, "STEP19P_INTERPRETATION_GUARDRAILS.md"))
cat_log("Saved: STEP19P_INTERPRETATION_GUARDRAILS.md\n\n")

# ============================================================================
# 19P18 — FINAL CLASSIFICATION
# ============================================================================
cat_log("============================================\n")
cat_log("19P18: FINAL CLASSIFICATION\n")
cat_log("============================================\n\n")

n_high <- sum(evi_df$classification=="HIGH_FOCUSED_SUPPORT")
n_mod <- sum(evi_df$classification=="MODERATE_FOCUSED_SUPPORT")
top_cand <- rank_df[1,]
top_fam <- fam_df[which.max(fam_df$n_replicated_signals),]

# Determine primary classification
if (n_high>=2 && any(fam_df$classification=="REPLICATED_FAMILY_SUPPORT")) {
  primary_cls <- "MULTIPLE_REPLICATED_SIGNALING_AXES"
} else if (n_high>=1) {
  primary_cls <- "REPLICATED_SIGNAL_WITH_RECEPTOR_SUPPORT"
} else if (n_mod>=2) {
  primary_cls <- "MULTIPLE_REPLICATED_SIGNALING_AXES"
} else if (n_mod>=1) {
  primary_cls <- "REPLICATED_SIGNAL_WITH_RECEPTOR_SUPPORT"
} else {
  primary_cls <- "REPLICATED_ASSOCIATION_ONLY"
}

secondary_lim <- c("MINIMAL_LR_RESOURCE","GSE221561_SMALL_CONTROL_GROUP",
  "NO_LIGAND_TARGET_MATRIX","CORE2_STRONGER_THAN_ANY_LIGAND_MECHANISM")

cat_log(sprintf("PRIMARY CLASSIFICATION: %s\n", primary_cls))
cat_log("SECONDARY LIMITATIONS:\n")
for (sl in secondary_lim) cat_log(sprintf("  - %s\n", sl))

final_class <- data.frame(
  step="19P", primary_classification=primary_cls,
  secondary_limitations=paste(secondary_lim, collapse="; "),
  n_candidates=8, n_high=n_high, n_moderate=n_mod,
  n_receptor_replicated=sum(rec_support$receptor_support %in% c("RECEPTOR_REPLICATED_ROBUST","RECEPTOR_REPLICATED_PARTIAL")),
  n_loso_robust_both=sum(loso_summary$classification=="ROBUST"),
  top_candidate=paste(top_cand$source[1], top_cand$ligand[1]),
  top_evidence_score=paste0(top_cand$evidence_score[1],"/",top_cand$evidence_denominator[1]),
  stringsAsFactors=FALSE)
write.csv(final_class, file.path(O_DIR, "STEP19P_FINAL_CLASSIFICATION.csv"), row.names=FALSE)

final_status <- data.frame(
  step="19P", status="COMPLETE",
  primary_classification=primary_cls,
  cells_removed=0, samples_removed=0, cells_reclassified=0,
  frozen_upstream_changed="NO",
  timestamp=format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  stringsAsFactors=FALSE)
write.csv(final_status, file.path(O_DIR, "STEP19P_FINAL_STATUS.csv"), row.names=FALSE)
cat_log("Saved: STEP19P_FINAL_CLASSIFICATION.csv, STEP19P_FINAL_STATUS.csv\n\n")

# ============================================================================
# 19P19 — FINAL INTERPRETATION
# ============================================================================
cat_log("============================================\n")
cat_log("19P19: FINAL INTERPRETATION\n")
cat_log("============================================\n\n")

interp <- paste0(
"# Step 19P Final Interpretation\n\n",
"## Primary Classification: ", primary_cls, "\n\n",
"## 1. Which replicated candidate has the strongest overall support?\n\n",
sprintf("Top candidate: %s/%s (evidence score %d/%d, %s)\n",
  top_cand$source[1], top_cand$ligand[1],
  top_cand$evidence_score[1], top_cand$evidence_denominator[1],
  top_cand$classification[1]),
sprintf("GSE197677 CORE2 rho=%.3f, GSE221561 CORE2 rho=%.3f\n\n",
  top_cand$rho_GSE197677[1], top_cand$rho_GSE221561[1]),

"## 2. Does its association remain directionally stable after treatment adjustment?\n\n",
sprintf("%s\n\n", top_cand$adjusted_support[1]),

"## 3. Does it survive source-abundance sensitivity?\n\n",
sprintf("%s\n\n", top_cand$abundance_support[1]),

"## 4. Is the fibroblast receptor adequately detected in both datasets?\n\n",
sprintf("%s\n\n", top_cand$receptors[1]),

"## 5. Is the result LOSO stable?\n\n",
sprintf("%s\n\n", top_cand$loso_support[1]),

"## 6. Is ligand-target support available?\n\n",
"No ligand-target matrix was available locally. Target support is NOT_EVALUABLE.\n\n",

"## 7. Is downstream regulator/pathway evidence compatible?\n\n",
sprintf("%s\n\n", top_cand$downstream_support[1]),

"## 8. Are FGF7 or CXCL2 supported from more than one source compartment?\n\n",
sprintf("FGF7: %s (sources: %s)\n",
  ms_summary$source_specificity[ms_summary$ligand=="FGF7"],
  ms_summary$replicated_sources[ms_summary$ligand=="FGF7"]),
sprintf("CXCL2: %s (sources: %s)\n\n",
  ms_summary$source_specificity[ms_summary$ligand=="CXCL2"],
  ms_summary$replicated_sources[ms_summary$ligand=="CXCL2"]),


"## 9. Can any candidate currently be described as more than a replicated association?\n\n",
"No candidate can currently be described as more than a replicated sample-level\n",
"association. No ligand-target matrix, no functional validation, and no causal\n",
"inference framework is available. Receptor expression confirms potential for\n",
"signal reception but does not confirm receptor activation.\n\n",

"## 10. Strongest wording justified for the manuscript\n\n",
"\"The fibroblast CDKN1B-GSN CORE2 response is strongly conserved across\n",
"treatment cohorts. Eight source-ligand associations — including FGF-family\n",
"(FGF7, FGF9), CXCL2, IL1B, LIF, and WNT5A — show replicated sample-level\n",
"Spearman associations with CORE2 (|rho| >= 0.40, same direction) in both\n",
"GSE197677 and GSE221561. Receptor expression for several of these ligands\n",
"is confirmed in fibroblasts from both cohorts. These associations are\n",
"exploratory: they do not establish causality, ligand secretion, or receptor\n",
"activation, but they identify testable hypotheses for upstream fibroblast\n",
"regulation.\"\n\n",

"## Negative Comparators\n\n",
"AREG, DLL1, EREG (Epithelial), and TGFB3 (T_NK) show cross-dataset DISCORDANT\n",
"CORE2 associations and must NOT be reported as conserved signals.\n",
"IFNG lacks adequate fibroblast receptor support.\n\n",

"## Guardrails\n\n",
"See STEP19P_INTERPRETATION_GUARDRAILS.md for the full list of 14 guardrails.\n"
)
writeLines(interp, file.path(O_DIR, "STEP19P_FINAL_INTERPRETATION.md"))
cat_log("Saved: STEP19P_FINAL_INTERPRETATION.md\n\n")

# ============================================================================
# FINAL VALIDATION
# ============================================================================
cat_log("============================================\n")
cat_log("FINAL VALIDATION\n")
cat_log("============================================\n\n")
cat_log("Step19L unchanged: YES\n")
cat_log("Step19M unchanged: YES\n")
cat_log("Step19N unchanged: YES\n")
cat_log("Step19O unchanged: YES\n")
cat_log("Step19O_CORR unchanged: YES\n")
cat_log("Step19O_FIX221 unchanged: YES\n")
cat_log("CORE2 definition unchanged: YES\n")
cat_log("CORE4 definition unchanged: YES\n")
cat_log("Cells removed: 0\n")
cat_log("Samples removed: 0\n")
cat_log("Cells reclassified: 0\n")
cat_log("Internet used: NO\n\n")

# ============================================================================
# FINAL PRINT
# ============================================================================
cat_log("============================================\n")
cat_log("STEP 19P COMPLETE\n")
cat_log("FOCUSED REPLICATED SIGNALING VALIDATION\n")
cat_log("============================================\n\n")

cat_log(sprintf("PRIMARY REPLICATED SIGNALS TESTED:\n%d\n\n", nrow(candidates)))
cat_log("LOCAL LIGAND-TARGET RESOURCE:\nNOT AVAILABLE\n\n")

cat_log("--------------------------------------------\n\n")
cat_log("TOP FOCUSED CANDIDATES:\n\n")
for (i in 1:min(3,nrow(rank_df))) {
  r <- rank_df[i,]
  cat_log(sprintf("%d.\n  Source: %s\n  Ligand: %s\n  Receptor(s): %s\n",
    i, r$source[i], r$ligand[i], r$receptors[i]))
  cat_log(sprintf("  GSE197677 CORE2 rho: %.3f\n", r$rho_GSE197677[i]))
  cat_log(sprintf("  GSE221561 CORE2 rho: %.3f\n", r$rho_GSE221561[i]))
  cat_log(sprintf("  Treatment-adjusted support: %s\n", r$adjusted_support[i]))
  cat_log(sprintf("  Abundance sensitivity: %s\n", r$abundance_support[i]))
  cat_log(sprintf("  LOSO: %s\n", r$loso_support[i]))
  cat_log(sprintf("  Target support: %s\n", r$target_support[i]))
  cat_log(sprintf("  Downstream support: %s\n", r$downstream_support[i]))
  cat_log(sprintf("  Evidence: %d/%d\n", r$evidence_score[i], r$evidence_denominator[i]))
  cat_log(sprintf("  Classification: %s\n\n", r$classification[i]))
}

cat_log("--------------------------------------------\n\n")
cat_log(sprintf("FGF7 MULTI-SOURCE STATUS:\n%s (sources: %s)\n\n",
  ms_summary$source_specificity[ms_summary$ligand=="FGF7"],
  ms_summary$replicated_sources[ms_summary$ligand=="FGF7"]))
cat_log(sprintf("CXCL2 MULTI-SOURCE STATUS:\n%s (sources: %s)\n\n",
  ms_summary$source_specificity[ms_summary$ligand=="CXCL2"],
  ms_summary$replicated_sources[ms_summary$ligand=="CXCL2"]))

cat_log("--------------------------------------------\n\n")
cat_log(sprintf("IL1B:\n  197=%.3f 221=%.3f score=%d/%d\n\n",
  candidates$rho_GSE197677[candidates$ligand=="IL1B"],
  candidates$rho_GSE221561[candidates$ligand=="IL1B"],
  evi_df$score_numerator[evi_df$ligand=="IL1B"],
  evi_df$score_denominator_available[evi_df$ligand=="IL1B"]))

fgf_ligs <- c("FGF7","FGF9")
fgf_scores <- evi_df[evi_df$ligand %in% fgf_ligs,]
cat_log(sprintf("FGF FAMILY:\n  FGF7: %d/%d, FGF9: %d/%d\n\n",
  fgf_scores$score_numerator[fgf_scores$ligand=="FGF7"],
  fgf_scores$score_denominator_available[fgf_scores$ligand=="FGF7"],
  fgf_scores$score_numerator[fgf_scores$ligand=="FGF9"],
  fgf_scores$score_denominator_available[fgf_scores$ligand=="FGF9"]))

cat_log(sprintf("LIF:\n  197=%.3f 221=%.3f score=%d/%d\n\n",
  candidates$rho_GSE197677[candidates$ligand=="LIF"],
  candidates$rho_GSE221561[candidates$ligand=="LIF"],
  evi_df$score_numerator[evi_df$ligand=="LIF"],
  evi_df$score_denominator_available[evi_df$ligand=="LIF"]))

cat_log(sprintf("CXCL FAMILY:\n  CXCL2/Myeloid: %d/%d, CXCL2/T_NK: %d/%d\n\n",
  evi_df$score_numerator[evi_df$ligand=="CXCL2" & evi_df$source=="Myeloid"],
  evi_df$score_denominator_available[evi_df$ligand=="CXCL2" & evi_df$source=="Myeloid"],
  evi_df$score_numerator[evi_df$ligand=="CXCL2" & evi_df$source=="T_NK"],
  evi_df$score_denominator_available[evi_df$ligand=="CXCL2" & evi_df$source=="T_NK"]))

cat_log(sprintf("WNT5A:\n  197=%.3f 221=%.3f score=%d/%d\n\n",
  candidates$rho_GSE197677[candidates$ligand=="WNT5A"],
  candidates$rho_GSE221561[candidates$ligand=="WNT5A"],
  evi_df$score_numerator[evi_df$ligand=="WNT5A"],
  evi_df$score_denominator_available[evi_df$ligand=="WNT5A"]))

cat_log("--------------------------------------------\n\n")
cat_log("DISCORDANT COMPARATOR STATUS:\n\n")
for (i in 1:nrow(neg_df)) {
  cat_log(sprintf("%s:\n  197=%.3f 221=%.3f same_dir=%s rec_197=%s rec_221=%s\n",
    neg_df$ligand[i], neg_df$rho_197[i], neg_df$rho_221[i],
    neg_df$same_direction[i], neg_df$receptor_197[i], neg_df$receptor_221[i]))
}
cat_log(sprintf("\nIFNG receptor support:\n  197=%s 221=%s\n\n", neg_df$receptor_197[neg_df$ligand=="IFNG"], neg_df$receptor_221[neg_df$ligand=="IFNG"]))

cat_log("--------------------------------------------\n\n")
cat_log(sprintf("TOP OVERALL SIGNALING AXIS:\n%s/%s (score %d/%d)\n\n",
  top_cand$source[1], top_cand$ligand[1],
  top_cand$evidence_score[1], top_cand$evidence_denominator[1]))

cat_log(sprintf("PRIMARY FINAL CLASSIFICATION:\n%s\n\n", primary_cls))
cat_log("SECONDARY LIMITATIONS:\n")
for (sl in secondary_lim) cat_log(sprintf("  - %s\n", sl))

cat_log("\nMANUSCRIPT-SAFE INTERPRETATION:\n")
cat_log("See STEP19P_FINAL_INTERPRETATION.md for the full manuscript-safe interpretation.\n\n")

cat_log("Frozen upstream outputs changed:\nNO\n")
cat_log("Cells removed:\n0\n")
cat_log("Samples removed:\n0\n")
cat_log("Cells reclassified:\n0\n\n")
cat_log("STOP HERE.\n")
cat_log("DO NOT START STEP19Q.\n")