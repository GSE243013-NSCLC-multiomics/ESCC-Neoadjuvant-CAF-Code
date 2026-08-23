# ============================================================
# ESCC Neoadjuvant Project
# STEP 13C
# GSE160269 metadata <-> barcode matching audit
#
# ONLY AUDIT:
# - no object modification
# - no merge
# - no QC
# - no normalization
# - no integration
# ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
})

cat("\n============================================================\n")
cat("STEP 13C: GSE160269 METADATA MATCH AUDIT\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------

raw_dir <- "01_raw/GEO/GSE160269"
obj_dir <- "03_objects/GSE160269/parts"
out_dir <- "06_tables/GSE160269/metadata"

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 1. Find 8 Seurat part objects
# ------------------------------------------------------------

object_files <- list.files(
  obj_dir,
  pattern = "\\.rds$",
  full.names = TRUE
)

if (length(object_files) != 8) {
  stop(
    "Expected 8 part objects, found ",
    length(object_files)
  )
}

cat("===== 1. READING CELL INDEX FROM 8 OBJECTS =====\n\n")

cell_index_list <- list()

for (f in sort(object_files)) {

  cat("Reading:", basename(f), "\n")

  obj <- readRDS(f)

  md <- obj[[]]

  if (!"original_barcode" %in% colnames(md)) {
    stop(
      "original_barcode missing in ",
      basename(f)
    )
  }

  if (!"source_celltype" %in% colnames(md)) {
    stop(
      "source_celltype missing in ",
      basename(f)
    )
  }

  tmp <- data.table(
    current_cell_id = colnames(obj),
    original_barcode = as.character(
      obj$original_barcode
    ),
    source_celltype = as.character(
      obj$source_celltype
    ),
    object_file = basename(f)
  )

  cell_index_list[[basename(f)]] <- tmp

  cat(
    "  cells:",
    nrow(tmp),
    "\n"
  )

  rm(obj, md)
  gc(verbose = FALSE)
}

cell_index <- rbindlist(
  cell_index_list,
  fill = TRUE
)

cat("\nTotal cells:", nrow(cell_index), "\n")
cat(
  "Expected   : 208659\n"
)

if (nrow(cell_index) != 208659) {
  stop(
    "Cell total does not equal 208659."
  )
}

cat("✓ Total cell count = 208659\n")

# ------------------------------------------------------------
# 2. Barcode duplication audit
# ------------------------------------------------------------

barcode_dup <- cell_index[
  ,
  .N,
  by = original_barcode
][
  N > 1
]

cat("\n===== 2. ORIGINAL BARCODE DUPLICATION =====\n\n")

cat(
  "Unique original barcodes:",
  uniqueN(cell_index$original_barcode),
  "\n"
)

cat(
  "Barcode strings occurring >1 time:",
  nrow(barcode_dup),
  "\n"
)

if (nrow(barcode_dup) > 0) {

  cat(
    "\nNOTE: repeated barcode strings are allowed across ",
    "different source matrices.\n"
  )

  fwrite(
    barcode_dup,
    file.path(
      out_dir,
      "GSE160269_repeated_original_barcodes.csv"
    ),
    bom = TRUE
  )
}

# ------------------------------------------------------------
# 3. Find probable author metadata files
# ------------------------------------------------------------

raw_files <- list.files(
  raw_dir,
  recursive = TRUE,
  full.names = TRUE
)

meta_files <- raw_files[
  grepl(
    "meta|annot",
    basename(raw_files),
    ignore.case = TRUE
  ) &
  grepl(
    "\\.(csv|tsv|txt)(\\.gz)?$",
    basename(raw_files),
    ignore.case = TRUE
  )
]

cat("\n===== 3. AUTHOR METADATA FILES =====\n\n")

if (length(meta_files) == 0) {
  stop(
    "No probable metadata files detected."
  )
}

for (f in meta_files) {

  cat(
    basename(f),
    " | ",
    round(file.info(f)$size / 1024^2, 2),
    " MB\n",
    sep = ""
  )
}

# ------------------------------------------------------------
# Helper: remove our cell-type prefix if one is present
# ------------------------------------------------------------

labels <- sort(
  unique(cell_index$source_celltype),
  decreasing = TRUE
)

strip_known_prefix <- function(x) {

  x <- trimws(as.character(x))

  for (lab in labels) {

    x <- sub(
      paste0(
        "^",
        lab,
        "(__|_|:|-)"
      ),
      "",
      x,
      ignore.case = TRUE
    )
  }

  x
}

# Additional normalization: strip trailing -1 suffix (10x style)
strip_10x_suffix <- function(x) {
  sub("-1$", "", x)
}

# Extract core barcode from sample-prefixed format
# e.g., P1T-I-AAACGGGAGCGAAGGG -> AAACGGGAGCGAAGGG
extract_core_barcode <- function(x) {
  # Pattern: SAMPLE-LETTERS-BARCODE
  sub("^[A-Z0-9]+-[A-Z]+-", "", x)
}

all_original <- unique(
  cell_index$original_barcode
)

all_current <- unique(
  cell_index$current_cell_id
)

# Also create set of core barcodes (prefix stripped)
all_core <- unique(
  extract_core_barcode(all_original)
)

current_to_original <- setNames(
  cell_index$original_barcode,
  cell_index$current_cell_id
)

# ------------------------------------------------------------
# 4. Read each metadata file and score EVERY column
# ------------------------------------------------------------

cat("\n===== 4. METADATA COLUMN SCORING =====\n\n")

metadata_objects <- list()
score_tables <- list()
selected_columns <- list()

for (f in meta_files) {

  cat(
    "\n------------------------------------------------------------\n"
  )

  cat(
    "Metadata file: ",
    basename(f),
    "\n",
    sep = ""
  )

  m <- tryCatch(
    fread(
      f,
      encoding = "UTF-8"
    ),
    error = function(e) {
      stop(
        "Could not read ",
        f,
        ": ",
        conditionMessage(e)
      )
    }
  )

  cat(
    "Rows:",
    nrow(m),
    " | Columns:",
    ncol(m),
    "\n"
  )

  scores <- rbindlist(
    lapply(
      names(m),
      function(col) {

        x <- trimws(
          as.character(
            m[[col]]
          )
        )

        stripped <- strip_known_prefix(x)
        no_10x <- strip_10x_suffix(x)
        core <- extract_core_barcode(x)
        stripped_10x <- strip_10x_suffix(stripped)

        data.table(
          metadata_file = basename(f),
          column = col,

          exact_original =
            sum(x %in% all_original, na.rm = TRUE),

          exact_current =
            sum(x %in% all_current, na.rm = TRUE),

          stripped_original =
            sum(stripped %in% all_original, na.rm = TRUE),

          suffix_10x_original =
            sum(no_10x %in% all_original, na.rm = TRUE),

          suffix_10x_core =
            sum(no_10x %in% all_core, na.rm = TRUE),

          core_original =
            sum(core %in% all_core, na.rm = TRUE),

          stripped_10x_original =
            sum(stripped_10x %in% all_original, na.rm = TRUE)
        )
      }
    )
  )

  scores[
    ,
    best_score :=
      pmax(
        exact_original,
        exact_current,
        stripped_original,
        suffix_10x_original,
        suffix_10x_core,
        core_original,
        stripped_10x_original
      )
  ]

  setorder(
    scores,
    -best_score
  )

  print(
    head(
      scores,
      10
    )
  )

  best_col <- scores$column[1]
  best_score <- scores$best_score[1]

  if (best_score == 0) {
    stop(
      "No barcode-like metadata column detected in ",
      basename(f)
    )
  }

  cat(
    "\nSelected barcode field: ",
    best_col,
    "\n",
    sep = ""
  )

  cat(
    "Best match score      : ",
    best_score,
    "\n",
    sep = ""
  )

  metadata_objects[[basename(f)]] <- m

  score_tables[[basename(f)]] <- scores

  selected_columns[[basename(f)]] <- data.table(
    metadata_file = basename(f),
    barcode_column = best_col,
    best_score = best_score
  )
}

score_all <- rbindlist(
  score_tables,
  fill = TRUE
)

selected_all <- rbindlist(
  selected_columns,
  fill = TRUE
)

fwrite(
  score_all,
  file.path(
    out_dir,
    "GSE160269_metadata_column_scores.csv"
  ),
  bom = TRUE
)

fwrite(
  selected_all,
  file.path(
    out_dir,
    "GSE160269_selected_barcode_fields.csv"
  ),
  bom = TRUE
)

# ------------------------------------------------------------
# 5. Build normalized metadata barcode pool
# ------------------------------------------------------------

cat("\n===== 5. BUILDING AUTHOR BARCODE POOL =====\n\n")

meta_barcode_list <- list()

for (nm in names(metadata_objects)) {

  m <- metadata_objects[[nm]]

  selected_col <- selected_all[
    metadata_file == nm,
    barcode_column
  ][1]

  raw_id <- trimws(
    as.character(
      m[[selected_col]]
    )
  )

  stripped_id <- strip_known_prefix(raw_id)

  # Additional normalizations for TCR-style barcodes
  # e.g., AAACCTGAGACAAAGG-1 -> AAACCTGAGACAAAGG
  no_10x_suffix <- strip_10x_suffix(raw_id)
  no_10x_suffix_stripped <- strip_10x_suffix(stripped_id)

  # Extract core barcode from sample-prefixed format
  # e.g., P1T-I-AAACGGGAGCGAAGGG -> AAACGGGAGCGAAGGG
  core_id <- extract_core_barcode(raw_id)

  normalized <- rep(
    NA_character_,
    length(raw_id)
  )

  mode <- rep(
    "unmatched",
    length(raw_id)
  )

  # 1. Exact original barcode
  hit <- raw_id %in% all_original
  normalized[hit] <- raw_id[hit]
  mode[hit] <- "exact_original"

  # 2. Exact current Seurat cell name
  hit2 <- (
    is.na(normalized) &
    raw_id %in% all_current
  )
  normalized[hit2] <- current_to_original[raw_id[hit2]]
  mode[hit2] <- "exact_current"

  # 3. Known cell-type prefix removed
  hit3 <- (
    is.na(normalized) &
    stripped_id %in% all_original
  )
  normalized[hit3] <- stripped_id[hit3]
  mode[hit3] <- "prefix_stripped"

  # 4. 10x suffix removed (e.g., -1)
  hit4 <- (
    is.na(normalized) &
    no_10x_suffix %in% all_original
  )
  normalized[hit4] <- no_10x_suffix[hit4]
  mode[hit4] <- "10x_suffix_stripped"

  # 5. Suffix 10x core: strip -1 THEN check all_core
  #    e.g., AAACCTGAGACAAAGG-1 -> AAACCTGAGACAAAGG -> check all_core
  no_10x_core <- strip_10x_suffix(raw_id)
  core_to_original <- setNames(
    all_original,
    extract_core_barcode(all_original)
  )
  hit5 <- (
    is.na(normalized) &
    no_10x_core %in% all_core
  )
  normalized[hit5] <- core_to_original[no_10x_core[hit5]]
  mode[hit5] <- "suffix_10x_core"

  # 6. Core barcode extracted (e.g., from P1T-I-XXX)
  hit6 <- (
    is.na(normalized) &
    core_id %in% all_core
  )
  normalized[hit6] <- core_to_original[core_id[hit6]]
  mode[hit6] <- "core_extracted"

  # 7. Both prefix stripped + 10x suffix stripped
  hit7 <- (
    is.na(normalized) &
    no_10x_suffix_stripped %in% all_original
  )
  normalized[hit7] <- no_10x_suffix_stripped[hit7]
  mode[hit7] <- "prefix_and_10x_stripped"

  tmp <- data.table(
    metadata_file = nm,
    metadata_row = seq_len(nrow(m)),
    raw_cell_id = raw_id,
    normalized_barcode = normalized,
    match_mode = mode
  )

  meta_barcode_list[[nm]] <- tmp

  cat(
    nm,
    ": matched ",
    sum(!is.na(normalized)),
    " / ",
    length(normalized),
    "\n",
    sep = ""
  )
}

meta_barcode_pool <- rbindlist(
  meta_barcode_list,
  fill = TRUE
)

fwrite(
  meta_barcode_pool,
  file.path(
    out_dir,
    "GSE160269_metadata_barcode_pool.csv.gz"
  )
)

# ------------------------------------------------------------
# 6. Match rate by each of the 8 Seurat objects
# ------------------------------------------------------------

cat("\n===== 6. BARCODE MATCH RATE BY CELL TYPE =====\n\n")

available_barcodes <- unique(
  meta_barcode_pool[
    !is.na(normalized_barcode),
    normalized_barcode
  ]
)

part_audit <- cell_index[
  ,
  .(
    matrix_cells = .N,

    barcode_found_in_author_metadata =
      sum(
        original_barcode %in%
          available_barcodes
      )
  ),
  by = source_celltype
]

part_audit[
  ,
  unmatched_cells :=
    matrix_cells -
    barcode_found_in_author_metadata
]

part_audit[
  ,
  match_percent :=
    round(
      100 *
      barcode_found_in_author_metadata /
      matrix_cells,
      3
    )
]

setorder(
  part_audit,
  source_celltype
)

print(part_audit)

fwrite(
  part_audit,
  file.path(
    out_dir,
    "GSE160269_barcode_match_by_celltype.csv"
  ),
  bom = TRUE
)

# ------------------------------------------------------------
# 7. Overall match
# ------------------------------------------------------------

total_cells <- sum(
  part_audit$matrix_cells
)

total_matched <- sum(
  part_audit$barcode_found_in_author_metadata
)

overall_match <- (
  total_matched /
  total_cells
) * 100

cat("\n============================================================\n")
cat("OVERALL MATCH AUDIT\n")
cat("============================================================\n\n")

cat(
  "Matrix cells : ",
  total_cells,
  "\n",
  sep = ""
)

cat(
  "Matched cells: ",
  total_matched,
  "\n",
  sep = ""
)

cat(
  "Unmatched    : ",
  total_cells - total_matched,
  "\n",
  sep = ""
)

cat(
  sprintf(
    "Match rate   : %.3f%%\n",
    overall_match
  )
)

# ------------------------------------------------------------
# 8. Per-file metadata summary
# ------------------------------------------------------------

file_summary <- meta_barcode_pool[
  ,
  .(
    metadata_rows = .N,
    matched_rows =
      sum(
        !is.na(normalized_barcode)
      ),
    unmatched_rows =
      sum(
        is.na(normalized_barcode)
      )
  ),
  by = metadata_file
]

file_summary[
  ,
  matched_percent :=
    round(
      100 *
      matched_rows /
      metadata_rows,
      3
    )
]

cat("\n===== 7. METADATA FILE SUMMARY =====\n\n")

print(file_summary)

fwrite(
  file_summary,
  file.path(
    out_dir,
    "GSE160269_metadata_file_summary.csv"
  ),
  bom = TRUE
)

# ------------------------------------------------------------
# 9. Save lightweight cell index
# ------------------------------------------------------------

fwrite(
  cell_index,
  file.path(
    out_dir,
    "GSE160269_Seurat_cell_index.csv.gz"
  )
)

# ------------------------------------------------------------
# FINAL
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STEP 13C FINISHED ✓\n")
cat("============================================================\n\n")

cat("This was an AUDIT ONLY.\n\n")

cat("No Seurat object modified ✓\n")
cat("No objects merged ✓\n")
cat("No cells filtered ✓\n")
cat("No normalization ✓\n")
cat("No integration ✓\n")
cat("STEP 14 NOT RUN ✓\n")

cat("\nSTOP HERE.\n")
