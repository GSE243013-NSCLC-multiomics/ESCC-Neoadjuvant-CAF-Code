# ============================================================
# ESCC Neoadjuvant Project
# STEP 17B
#
# CNV / malignant epithelial PRE-FLIGHT
#
# Checks:
#   1. candidate R packages
#   2. raw RNA counts availability
#   3. Seurat v5 count-layer structure
#   4. gene ID format
#   5. hg38 overlap if infercna is already installed
#   6. epithelial / normal-reference sample capacity
#   7. approximate memory requirements
#
# IMPORTANT:
#   NO package installation
#   NO JoinLayers
#   NO normalization
#   NO CNV inference
#   NO malignant calling
#   NO cell filtering
# ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
})

cat("\n============================================================\n")
cat("STEP 17B: CNV PREFLIGHT\n")
cat("============================================================\n\n")

outdir <- "06_tables/cross_dataset/malignant_epithelial"

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# 1. SYSTEM / R ENVIRONMENT
# ============================================================

cat("===== 1. SYSTEM / R ENVIRONMENT =====\n\n")

cat(
  "R version      : ",
  R.version.string,
  "\n",
  sep = ""
)

cat(
  "Platform       : ",
  R.version$platform,
  "\n",
  sep = ""
)

architecture <- tryCatch(
  system(
    "uname -m",
    intern = TRUE
  ),
  error = function(e) NA_character_
)

cat(
  "Architecture   : ",
  architecture,
  "\n",
  sep = ""
)

ram_bytes <- tryCatch(
  as.numeric(
    system(
      "sysctl -n hw.memsize",
      intern = TRUE
    )
  ),
  error = function(e) NA_real_
)

ram_gb <- ram_bytes / 1024^3

cat(
  "Physical RAM   : ",
  ifelse(
    is.finite(ram_gb),
    sprintf("%.1f GB", ram_gb),
    "unknown"
  ),
  "\n",
  sep = ""
)

cat("\nDisk:\n")

try(
  system(
    "df -h ."
  )
)

# ============================================================
# 2. PACKAGE STATUS
# ============================================================

cat("\n===== 2. CNV PACKAGE STATUS =====\n\n")

packages <- c(
  "infercna",
  "copykat",
  "infercnv",
  "remotes",
  "devtools",
  "Seurat",
  "SeuratObject",
  "Matrix",
  "data.table"
)

pkg_status <- rbindlist(
  lapply(
    packages,
    function(pkg) {

      installed <- requireNamespace(
        pkg,
        quietly = TRUE
      )

      version <- if (installed) {
        as.character(
          packageVersion(pkg)
        )
      } else {
        NA_character_
      }

      data.table(
        package = pkg,
        installed = installed,
        version = version
      )
    }
  )
)

print(pkg_status)

fwrite(
  pkg_status,
  file.path(
    outdir,
    "STEP17B_package_status.csv"
  ),
  bom = TRUE
)

cat(
  "\nNOTE:\n",
  "This step does NOT install any CNV package.\n",
  sep = ""
)

# ============================================================
# 3. STEP 17A EPITHELIAL INVENTORY
# ============================================================

inventory_file <- file.path(
  outdir,
  "STEP17A_all_sample_epithelial_inventory.csv"
)

if (!file.exists(inventory_file)) {
  stop(
    "Missing Step 17A inventory: ",
    inventory_file
  )
}

epi_inventory <- fread(
  inventory_file
)

cat("\n===== 3. STEP 17A EPITHELIAL INVENTORY =====\n\n")

epi_summary <- epi_inventory[
  ,
  .(
    samples = .N,

    epithelial_cells =
      sum(
        epithelial_cells,
        na.rm = TRUE
      ),

    samples_ge20 =
      sum(
        epithelial_cells >= 20
      ),

    samples_ge30 =
      sum(
        epithelial_cells >= 30
      ),

    samples_ge100 =
      sum(
        epithelial_cells >= 100
      )
  ),
  by = .(
    dataset,
    tissue
  )
]

print(epi_summary)

normal_details <- epi_inventory[
  tissue == "Normal",
  .(
    dataset,
    patient_id,
    sample_id,
    epithelial_cells,
    ge20 = epithelial_cells >= 20,
    ge30 = epithelial_cells >= 30,
    ge100 = epithelial_cells >= 100
  )
]

cat("\nNormal epithelial samples:\n\n")

print(
  normal_details,
  nrows = 100
)

fwrite(
  normal_details,
  file.path(
    outdir,
    "STEP17B_normal_reference_sample_details.csv"
  ),
  bom = TRUE
)

# ============================================================
# 4. LOCATE THE THREE SOURCE OBJECTS
# ============================================================

cat("\n===== 4. LOCATING SOURCE OBJECTS =====\n\n")

find_object <- function(
  directory,
  preferred_patterns
) {

  files <- list.files(
    directory,
    pattern = "\\.rds$",
    full.names = TRUE,
    recursive = TRUE,
    ignore.case = TRUE
  )

  if (length(files) == 0) {
    return(NA_character_)
  }

  for (pat in preferred_patterns) {

    hit <- files[
      grepl(
        pat,
        basename(files),
        ignore.case = TRUE
      )
    ]

    if (length(hit) > 0) {
      return(hit[1])
    }
  }

  files[1]
}

objects <- data.table(

  dataset = c(
    "GSE221561",
    "GSE197677",
    "GSE160269"
  ),

  object_file = c(

    find_object(
      "03_objects/GSE221561",
      c(
        "author_annotated",
        "counts_merged"
      )
    ),

    find_object(
      "03_objects/GSE197677",
      c(
        "^GSE197677_integrated\\.rds$",
        "integrated"
      )
    ),

    find_object(
      "03_objects/GSE160269",
      c(
        "annotated_merged_raw",
        "merged"
      )
    )
  )
)

objects[
  ,
  file_exists :=
    file.exists(object_file)
]

objects[
  ,
  file_size_GB :=
    fifelse(
      file_exists,
      file.info(object_file)$size / 1024^3,
      NA_real_
    )
]

print(objects)

if (any(!objects$file_exists)) {

  stop(
    "At least one required Seurat object is missing."
  )
}

fwrite(
  objects,
  file.path(
    outdir,
    "STEP17B_object_inventory.csv"
  ),
  bom = TRUE
)

# ============================================================
# 5. infercna hg38 GENOME REFERENCE
# ============================================================

infercna_installed <- requireNamespace(
  "infercna",
  quietly = TRUE
)

hg38_symbols <- NULL

if (infercna_installed) {

  cat("\n===== 5. infercna hg38 GENOME =====\n\n")

  hg38 <- tryCatch(
    infercna::retrieveGenome(
      name = "hg38"
    ),
    error = function(e) {
      cat(
        "Could not retrieve infercna hg38 genome:\n",
        conditionMessage(e),
        "\n"
      )
      NULL
    }
  )

  if (!is.null(hg38)) {

    hg38_symbols <- unique(
      as.character(
        hg38$symbol
      )
    )

    cat(
      "hg38 genome symbols:",
      length(hg38_symbols),
      "\n"
    )
  }

} else {

  cat("\n===== 5. infercna hg38 GENOME =====\n\n")

  cat(
    "infercna is NOT installed.\n",
    "Gene-overlap calculation will be deferred.\n",
    sep = ""
  )
}

# ============================================================
# 6. FUNCTIONS FOR OBJECT AUDIT
# ============================================================

get_layers_safe <- function(assay_obj) {

  # Direct slot inspection — no SeuratObject::Layers call
  tryCatch({

    slot_names <- slotNames(assay_obj)

    # Look for names/counts slots
    count_slots <- slot_names[
      grepl(
        "^(counts|data|scale\\.data|names)$",
        slot_names,
        ignore.case = TRUE
      )
    ]

    # If there is a "names" slot, it may contain layer names
    if ("names" %in% slot_names) {

      nm <- tryCatch(
        slot(assay_obj, "names"),
        error = function(e) character(0)
      )

      if (length(nm) > 0) {
        return(as.character(nm))
      }
    }

    if (length(count_slots) > 0) {
      return(count_slots)
    }

    character(0)

  }, error = function(e) {
    character(0)
  })
}

# Safe wrapper for GetAssayData / LayerData
get_counts_safe <- function(obj, assay_name) {

  # Try v5 LayerData approach first
  tryCatch({

    rna <- obj[[assay_name]]

    layers <- get_layers_safe(rna)

    count_layers <- layers[
      grepl(
        "^counts",
        layers,
        ignore.case = TRUE
      )
    ]

    if (length(count_layers) > 0) {

      result <- tryCatch(
        SeuratObject::LayerData(
          rna,
          layer = count_layers[1]
        ),
        error = function(e) NULL
      )

      if (!is.null(result)) {
        return(list(
          accessible = TRUE,
          status = paste0(
            "LAYER_",
            count_layers[1]
          ),
          data = result
        ))
      }
    }

    NULL

  }, error = function(e) NULL)

  # Fallback: v4 GetAssayData
  tryCatch({

    result <- SeuratObject::GetAssayData(
      obj,
      assay = assay_name,
      slot = "counts"
    )

    if (!is.null(result) &&
        nrow(result) > 0 &&
        ncol(result) > 0) {

      return(list(
        accessible = TRUE,
        status = "V4_COUNTS_SLOT",
        data = result
      ))
    }

    NULL

  }, error = function(e) NULL)
}

audit_one_object <- function(
  dataset,
  object_file
) {

  tryCatch({

  cat("\n")
  cat("============================================================\n")
  cat("AUDITING ", dataset, "\n", sep = "")
  cat("============================================================\n\n")

  cat(
    "Loading: ",
    object_file,
    "\n",
    sep = ""
  )

  obj <- readRDS(
    object_file
  )

  md <- obj[[]]

  n_cells <- nrow(md)

  # Defensive assay name extraction
  assays <- character(0)

  tryCatch({

    assay_slot <- obj@assays

    if (is(assay_slot, "AssayMap")) {
      assays <- names(assay_slot)
    } else if (is(assay_slot, "list")) {
      assays <- names(assay_slot)
    } else {
      assays <- tryCatch(
        SeuratObject::AssayNames(obj),
        error = function(e) character(0)
      )
    }

  }, error = function(e) {
    assays <<- character(0)
  })

  cat(
    "Metadata cells : ",
    n_cells,
    "\n",
    sep = ""
  )

  cat(
    "Assays         : ",
    paste(
      assays,
      collapse = ", "
    ),
    "\n",
    sep = ""
  )

  # Try to access RNA assay directly from slot
  rna <- NULL

  tryCatch({

    rna <- obj@assays[["RNA"]]

  }, error = function(e) {

    cat(
      "  Direct slot access failed, trying obj[[RNA]]...\n"
    )

    tryCatch({
      rna <<- obj[["RNA"]]
    }, error = function(e2) {
      cat(
        "  obj[[RNA]] also failed:",
        conditionMessage(e2),
        "\n"
      )
    })
  })

  if (is.null(rna)) {

    result <- data.table(
      dataset = dataset,
      object_file = object_file,
      cells = n_cells,
      genes = NA_integer_,
      RNA_assay = FALSE,
      count_layer_status = "NO_RNA_ASSAY",
      count_layers = NA_character_,
      raw_counts_accessible = FALSE,
      ensembl_like_pct = NA_real_,
      duplicated_gene_names = NA_integer_,
      hg38_symbol_overlap_n = NA_integer_,
      hg38_symbol_overlap_pct = NA_real_,
      sparse_counts_size_GB = NA_real_
    )

    rm(obj, md)
    gc(verbose = FALSE)

    return(result)
  }

  genes <- tryCatch(
    rownames(rna),
    error = function(e) {
      tryCatch(
        slot(rna, "dimnames")[[1]],
        error = function(e2) character(0)
      )
    }
  )

  n_genes <- length(
    genes
  )

  # Use safe wrapper to detect count layers and accessibility
  counts_info <- tryCatch(
    get_counts_safe(obj, "RNA"),
    error = function(e) NULL
  )

  if (!is.null(counts_info) && counts_info$accessible) {

    raw_counts_accessible <- TRUE
    count_status <- counts_info$status
    count_layers <- counts_info$status

    sparse_size <- as.numeric(
      object.size(counts_info$data)
    ) / 1024^3

    rm(counts_info)

  } else {

    raw_counts_accessible <- FALSE
    count_status <- "COUNTS_NOT_ACCESSIBLE"
    count_layers <- NA_character_
    sparse_size <- NA_real_
  }

  ensembl_like <- grepl(
    "^ENSG[0-9]+",
    genes,
    ignore.case = TRUE
  )

  ensembl_pct <- 100 *
    mean(
      ensembl_like
    )

  duplicated_genes <- sum(
    duplicated(
      genes
    )
  )

  if (!is.null(hg38_symbols)) {

    overlap_n <- sum(
      unique(genes) %in%
        hg38_symbols
    )

    overlap_pct <- 100 *
      overlap_n /
      length(
        unique(genes)
      )

  } else {

    overlap_n <- NA_integer_
    overlap_pct <- NA_real_
  }

  cat(
    "RNA genes      : ",
    n_genes,
    "\n",
    sep = ""
  )

  cat(
    "Count status   : ",
    count_status,
    "\n",
    sep = ""
  )

  cat(
    "Count layers   : ",
    ifelse(
      is.na(count_layers),
      "<not detected>",
      count_layers
    ),
    "\n",
    sep = ""
  )

  cat(
    "ENSG-like genes: ",
    sprintf(
      "%.2f%%",
      ensembl_pct
    ),
    "\n",
    sep = ""
  )

  cat(
    "Duplicated gene names: ",
    duplicated_genes,
    "\n",
    sep = ""
  )

  if (!is.na(overlap_pct)) {

    cat(
      "hg38 symbol overlap: ",
      overlap_n,
      " / ",
      n_genes,
      " (",
      sprintf("%.1f%%", overlap_pct),
      ")\n",
      sep = ""
    )
  }

  result <- data.table(
    dataset = dataset,

    object_file =
      object_file,

    object_size_GB =
      file.info(object_file)$size /
      1024^3,

    cells =
      n_cells,

    genes =
      n_genes,

    RNA_assay =
      TRUE,

    count_layer_status =
      count_status,

    count_layers =
      ifelse(
        length(count_layers) == 0,
        NA_character_,
        paste(
          count_layers,
          collapse = "; "
        )
      ),

    raw_counts_accessible =
      raw_counts_accessible,

    ensembl_like_pct =
      round(
        ensembl_pct,
        3
      ),

    duplicated_gene_names =
      duplicated_genes,

    hg38_symbol_overlap_n =
      overlap_n,

    hg38_symbol_overlap_pct =
      round(
        overlap_pct,
        3
      ),

    sparse_counts_size_GB =
      round(
        sparse_size,
        3
      )
  )

  rm(
    obj,
    md,
    rna,
    genes
  )

  gc(verbose = FALSE)

  cat(
    "\n✓ ",
    dataset,
    " audit finished; memory released.\n",
    sep = ""
  )

  result

  }, error = function(e) {

    cat(
      "\n⚠ AUDIT FAILED for ",
      dataset,
      ": ",
      conditionMessage(e),
      "\n",
      sep = ""
    )

    data.table(
      dataset = dataset,
      object_file = object_file,
      cells = NA_integer_,
      genes = NA_integer_,
      RNA_assay = NA,
      count_layer_status = paste0(
        "AUDIT_ERROR: ",
        conditionMessage(e)
      ),
      count_layers = NA_character_,
      raw_counts_accessible = NA,
      ensembl_like_pct = NA_real_,
      duplicated_gene_names = NA_integer_,
      hg38_symbol_overlap_n = NA_integer_,
      hg38_symbol_overlap_pct = NA_real_,
      sparse_counts_size_GB = NA_real_
    )
  })
}

# ============================================================
# 7. AUDIT OBJECTS SEQUENTIALLY
# ============================================================

cat("\n===== 6. SEQUENTIAL OBJECT AUDIT =====\n")

audit_list <- list()

for (i in seq_len(nrow(objects))) {

  audit_list[[i]] <- audit_one_object(
    objects$dataset[i],
    objects$object_file[i]
  )
}

object_audit <- rbindlist(
  audit_list,
  fill = TRUE
)

cat("\n============================================================\n")
cat("7. OBJECT / GENE / COUNT AUDIT SUMMARY\n")
cat("============================================================\n\n")

print(
  object_audit
)

fwrite(
  object_audit,
  file.path(
    outdir,
    "STEP17B_object_gene_count_audit.csv"
  ),
  bom = TRUE
)

# ============================================================
# 8. EPITHELIAL MEMORY PLANNING
# ============================================================

cat("\n===== 8. EPITHELIAL MEMORY PLANNING =====\n\n")

epi_total <- epi_inventory[
  ,
  .(
    epithelial_cells =
      sum(
        epithelial_cells,
        na.rm = TRUE
      ),

    tumor_epi_cells =
      sum(
        epithelial_cells[
          tissue == "Tumor"
        ],
        na.rm = TRUE
      ),

    normal_epi_cells =
      sum(
        epithelial_cells[
          tissue == "Normal"
        ],
        na.rm = TRUE
      )
  ),
  by = dataset
]

memory_plan <- merge(
  object_audit[
    ,
    .(
      dataset,
      genes
    )
  ],
  epi_total,
  by = "dataset",
  all.x = TRUE
)

memory_plan[
  ,
  dense_epithelial_matrix_GB :=
    genes *
    epithelial_cells *
    8 /
    1024^3
]

memory_plan[
  ,
  rough_3x_working_GB :=
    dense_epithelial_matrix_GB * 3
]

memory_plan[
  ,
  physical_RAM_GB :=
    ram_gb
]

memory_plan[
  ,
  full_dataset_memory_flag :=
    fcase(

      !is.finite(
        physical_RAM_GB
      ),
      "RAM_UNKNOWN",

      rough_3x_working_GB <=
        physical_RAM_GB * 0.60,
      "LIKELY_MANAGEABLE",

      rough_3x_working_GB <=
        physical_RAM_GB * 0.90,
      "BORDERLINE_USE_SAMPLEWISE",

      default =
        "HIGH_RISK_USE_SAMPLEWISE"
    )
]

print(memory_plan)

cat(
  "\nNOTE:\n",
  "These are rough dense-matrix planning estimates only.\n",
  "They are intentionally conservative and are not runtime guarantees.\n",
  sep = ""
)

fwrite(
  memory_plan,
  file.path(
    outdir,
    "STEP17B_memory_planning.csv"
  ),
  bom = TRUE
)

# ============================================================
# 9. NORMAL REFERENCE GROUP FEASIBILITY
# ============================================================

cat("\n===== 9. NORMAL REFERENCE GROUP FEASIBILITY =====\n\n")

reference_plan <- epi_inventory[
  tissue == "Normal",
  .(
    normal_samples_total = .N,

    normal_samples_with_epi =
      sum(
        epithelial_cells > 0
      ),

    normal_samples_ge20 =
      sum(
        epithelial_cells >= 20
      ),

    normal_samples_ge30 =
      sum(
        epithelial_cells >= 30
      ),

    total_normal_epi =
      sum(
        epithelial_cells
      ),

    median_normal_epi =
      median(
        epithelial_cells
      )
  ),
  by = dataset
]

reference_plan[
  ,
  project_reference_flag :=
    fcase(

      normal_samples_ge30 >= 2,
      "STRONG_MULTI_SAMPLE_REFERENCE",

      normal_samples_ge20 >= 2,
      "USABLE_BUT_SMALL_REFERENCE",

      normal_samples_with_epi >= 2,
      "LIMITED_REFERENCE_SENSITIVITY_REQUIRED",

      default =
        "INSUFFICIENT_MULTI_SAMPLE_REFERENCE"
    )
]

print(reference_plan)

fwrite(
  reference_plan,
  file.path(
    outdir,
    "STEP17B_normal_reference_preflight.csv"
  ),
  bom = TRUE
)

# ============================================================
# 10. DATASET-SPECIFIC PREFLIGHT DECISION
# ============================================================

cat("\n============================================================\n")
cat("10. DATASET-SPECIFIC CNV PREFLIGHT DECISION\n")
cat("============================================================\n\n")

decision <- merge(
  object_audit,
  reference_plan,
  by = "dataset",
  all.x = TRUE
)

decision <- merge(
  decision,
  memory_plan[
    ,
    .(
      dataset,
      epithelial_cells,
      tumor_epi_cells,
      normal_epi_cells,
      dense_epithelial_matrix_GB,
      rough_3x_working_GB,
      full_dataset_memory_flag
    )
  ],
  by = "dataset",
  all.x = TRUE
)

decision[
  ,
  gene_format_flag :=
    fcase(

      ensembl_like_pct >= 50,
      "ENSEMBL_DOMINANT_NEEDS_MAPPING_CHECK",

      duplicated_gene_names > 0,
      "SYMBOLS_WITH_DUPLICATES_REVIEW",

      default =
        "GENE_SYMBOL_FORMAT_COMPATIBLE"
    )
]

decision[
  ,
  count_preparation_flag :=
    fcase(

      count_layer_status ==
        "MULTIPLE_COUNTS_LAYERS_NEED_JOIN_OR_SAMPLE_EXTRACTION",
      "DO_NOT_JOIN_GLOBAL_OBJECT_USE_SAMPLEWISE_EXTRACTION",

      raw_counts_accessible == TRUE,
      "RAW_COUNTS_AVAILABLE",

      default =
        "RAW_COUNTS_PROBLEM"
    )
]

decision[
  ,
  recommended_execution_mode :=
    fcase(

      dataset == "GSE160269",
      paste0(
        "SAMPLEWISE_CNV; pooled/balanced normal reference; ",
        "independent sensitivity analysis recommended"
      ),

      full_dataset_memory_flag %in%
        c(
          "HIGH_RISK_USE_SAMPLEWISE",
          "BORDERLINE_USE_SAMPLEWISE"
        ),
      "SAMPLEWISE_CNV_WITH_FIXED_REFERENCE",

      default =
        "SAMPLEWISE_CNV_WITH_FIXED_REFERENCE"
    )
]

decision[
  ,
  preflight_status :=
    fcase(

      raw_counts_accessible == FALSE,
      "BLOCKED_RAW_COUNTS",

      gene_format_flag ==
        "ENSEMBL_DOMINANT_NEEDS_MAPPING_CHECK",
      "REVIEW_GENE_IDS",

      project_reference_flag ==
        "INSUFFICIENT_MULTI_SAMPLE_REFERENCE",
      "REFERENCE_LIMITED",

      default =
        "READY_FOR_SMALL_PILOT"
    )
]

display_cols <- c(
  "dataset",
  "genes",
  "epithelial_cells",
  "tumor_epi_cells",
  "normal_epi_cells",
  "count_layer_status",
  "gene_format_flag",
  "hg38_symbol_overlap_pct",
  "project_reference_flag",
  "rough_3x_working_GB",
  "full_dataset_memory_flag",
  "recommended_execution_mode",
  "preflight_status"
)

print(
  decision[
    ,
    ..display_cols
  ],
  nrows = 100
)

fwrite(
  decision,
  file.path(
    outdir,
    "STEP17B_CNV_preflight_decision.csv"
  ),
  bom = TRUE
)

# ============================================================
# 11. METHOD STATUS
# ============================================================

cat("\n============================================================\n")
cat("11. METHOD STATUS\n")
cat("============================================================\n\n")

cat(
  "infercna installed : ",
  infercna_installed,
  "\n",
  sep = ""
)

cat(
  "copykat installed  : ",
  requireNamespace(
    "copykat",
    quietly = TRUE
  ),
  "\n",
  sep = ""
)

cat(
  "infercnv installed : ",
  requireNamespace(
    "infercnv",
    quietly = TRUE
  ),
  "\n",
  sep = ""
)

cat(
  "\nNo method is executed in Step 17B.\n"
)

# ============================================================
# 12. PROPOSED PILOT ORDER
# ============================================================

pilot_order <- data.table(

  order = 1:3,

  dataset = c(
    "GSE197677",
    "GSE221561",
    "GSE160269"
  ),

  rationale = c(
    paste0(
      "Smallest epithelial compartment (4,857 cells) ",
      "and 4 normal epithelial samples; safest pilot"
    ),

    paste0(
      "26,397 epithelial cells with 2 strong ",
      "adjacent-normal references"
    ),

    paste0(
      "Largest epithelial compartment (44,730 cells) ",
      "and limited normal-reference distribution; run last"
    )
  )
)

cat("\n============================================================\n")
cat("12. PROPOSED PILOT ORDER\n")
cat("============================================================\n\n")

print(pilot_order)

fwrite(
  pilot_order,
  file.path(
    outdir,
    "STEP17B_proposed_CNV_pilot_order.csv"
  ),
  bom = TRUE
)

# ============================================================
# FINAL
# ============================================================

cat("\n============================================================\n")
cat("STEP 17B FINISHED ✓\n")
cat("CNV PREFLIGHT ONLY ✓\n")
cat("============================================================\n\n")

cat("No packages installed ✓\n")
cat("No JoinLayers performed ✓\n")
cat("No expression normalization performed ✓\n")
cat("No CNV inference performed ✓\n")
cat("No malignant labels assigned ✓\n")
cat("No cells filtered ✓\n")
cat("No pseudobulk DE ✓\n")
cat("No cross-dataset integration ✓\n")

cat(
  "\nIMPORTANT:\n",
  "Do NOT start all three datasets at once.\n",
  "The next step should be ONE small pilot only.\n",
  sep = ""
)

cat("\nSTOP HERE.\n")
