#!/bin/bash

set -euo pipefail

PROJECT="$HOME/Downloads/ESCC Neoadjuvant"

RAW_DIR="$PROJECT/01_raw/GEO/GSE160269"
PROC_DIR="$PROJECT/02_processed/GSE160269/10x_parts"
OBJ_DIR="$PROJECT/03_objects/GSE160269/parts"
META_DIR="$PROJECT/00_metadata"
LOG_DIR="$PROJECT/08_logs"

mkdir -p "$PROC_DIR"
mkdir -p "$OBJ_DIR"
mkdir -p "$META_DIR"
mkdir -p "$LOG_DIR"

cd "$PROJECT"

echo "============================================================"
echo "STEP 13B: GSE160269 LOW-MEMORY BUILD"
echo ""
echo "ONLY:"
echo "  wide UMI matrix -> sparse MatrixMarket"
echo "  -> one Seurat object per cell type"
echo ""
echo "NO MERGE"
echo "NO QC"
echo "NO NORMALIZATION"
echo "NO INTEGRATION"
echo "NO STEP 14"
echo "============================================================"
echo ""


# ============================================================
# 1. Discover the 8 ordinary UMI matrices
#    Exclude CD45 / TCR / metadata
# ============================================================

echo "===== 1. Detecting GSE160269 UMI matrices ====="

python3 - <<'PY'
import os
import re
import sys

root = "01_raw/GEO/GSE160269"
out  = "00_metadata/06_GSE160269_UMI_matrix_manifest.tsv"

candidates = []

for dp, _, files in os.walk(root):

    for fn in files:

        low = fn.lower()

        # Must be a UMI matrix
        if "umi" not in low:
            continue

        # Do NOT include special CD45 matrices
        if "cd45" in low:
            continue

        # Exclude metadata/TCR
        if any(
            x in low
            for x in [
                "metadata",
                "meta",
                "tcr"
            ]
        ):
            continue

        if not (
            low.endswith(".txt") or
            low.endswith(".tsv") or
            low.endswith(".csv") or
            low.endswith(".txt.gz") or
            low.endswith(".tsv.gz") or
            low.endswith(".csv.gz")
        ):
            continue

        candidates.append(
            os.path.join(dp, fn)
        )

candidates = sorted(candidates)

print()
print("Detected candidate UMI matrices:")
print()

for i, p in enumerate(candidates, 1):
    size = os.path.getsize(p) / 1024**3
    print(
        f"{i:02d}. {os.path.basename(p)}"
        f"   ({size:.2f} GB)"
    )

print()
print("Number detected:", len(candidates))

if len(candidates) != 8:

    print()
    print("ERROR:")
    print("Expected exactly 8 non-CD45 UMI matrices.")
    print("Nothing will be built.")
    print()
    sys.exit(2)


def make_label(path):

    name = os.path.basename(path)

    # remove compression
    name = re.sub(
        r"\.gz$",
        "",
        name,
        flags=re.I
    )

    # remove text extension
    name = re.sub(
        r"\.(txt|tsv|csv)$",
        "",
        name,
        flags=re.I
    )

    # Extract cell type from pattern: GSE160269_UMI_matrix_CellType
    # Match the cell type portion after _UMI_matrix_
    match = re.search(
        r"_UMI_matrix_([^_]+)",
        name,
        flags=re.I
    )

    if match:
        return match.group(1)

    # Fallback: remove prefix and UMI keywords
    name = re.sub(
        r"^GSE160269[_\.-]*",
        "",
        name,
        flags=re.I
    )

    name = re.sub(
        r"(?i)_*(umi|matrix|_)*_*",
        "_",
        name
    )

    name = re.sub(
        r"[^A-Za-z0-9]+",
        "_",
        name
    )

    name = name.strip("_")

    return name or "unknown"


rows = []

labels_seen = set()

for p in candidates:

    label = make_label(p)

    if label in labels_seen:
        print(
            "ERROR: duplicate generated label:",
            label
        )
        sys.exit(3)

    labels_seen.add(label)

    rows.append(
        (label, os.path.abspath(p))
    )


with open(
    out,
    "w",
    encoding="utf-8"
) as f:

    f.write("label\tsource_file\n")

    for label, path in rows:
        f.write(
            f"{label}\t{path}\n"
        )


print()
print("Generated labels:")

for label, path in rows:
    print(
        f"  {label:25s} <- {os.path.basename(path)}"
    )

print()
print("Saved:", out)
PY


echo ""
echo "===== 2. Streaming conversion to sparse matrices ====="
echo ""
echo "This is intentionally memory-efficient."
echo "Each UMI file is processed line-by-line."
echo ""


# ============================================================
# 2. Streaming wide matrix -> MatrixMarket
# ============================================================

tail -n +2 \
"$META_DIR/06_GSE160269_UMI_matrix_manifest.tsv" |
while IFS=$'\t' read -r LABEL SOURCE
do

    OUT="$PROC_DIR/$LABEL"
    DONE="$OUT/.conversion_complete"

    mkdir -p "$OUT"

    echo ""
    echo "------------------------------------------------------------"
    echo "CELL TYPE / MATRIX: $LABEL"
    echo "SOURCE:"
    echo "$SOURCE"
    echo "------------------------------------------------------------"

    if [ -f "$DONE" ]; then
        echo "✓ Sparse conversion already completed — skipping."
        continue
    fi


    python3 - "$SOURCE" "$OUT" "$LABEL" <<'PY'
import sys
import os
import csv
import gzip
import shutil
import time

csv.field_size_limit(sys.maxsize)

source = sys.argv[1]
outdir = sys.argv[2]
label  = sys.argv[3]

os.makedirs(
    outdir,
    exist_ok=True
)

matrix_out = os.path.join(
    outdir,
    "matrix.mtx.gz"
)

feature_out = os.path.join(
    outdir,
    "features.tsv.gz"
)

barcode_out = os.path.join(
    outdir,
    "barcodes.tsv.gz"
)

tmp_entries = os.path.join(
    outdir,
    "matrix_entries.tmp"
)


def open_text(path):

    if path.lower().endswith(".gz"):
        return gzip.open(
            path,
            "rt",
            encoding="utf-8",
            errors="replace",
            newline=""
        )

    return open(
        path,
        "r",
        encoding="utf-8",
        errors="replace",
        newline=""
    )


print()
print("Opening:", os.path.basename(source))

# ------------------------------------------------------------
# Determine delimiter and header format
# GSE160269 files are SPACE-delimited:
#   Line 1: barcode1 barcode2 barcode3 ...
#   Line 2+: gene_name val1 val2 val3 ...
# ------------------------------------------------------------

with open_text(source) as fh:

    first_line = fh.readline()

    if not first_line:
        raise RuntimeError(
            "Empty UMI matrix"
        )

    tabs = first_line.count("\t")
    commas = first_line.count(",")
    spaces = len(first_line.split()) - 1

    # Detect delimiter: tab > comma > space
    if tabs > 0 and tabs >= commas:
        delimiter = "\t"
    elif commas > 0:
        delimiter = ","
    else:
        delimiter = " "  # space-delimited

    print(
        "Delimiter:",
        repr(delimiter)
    )

    if delimiter == " ":
        # Space-delimited: split on whitespace
        header = first_line.strip().split()
    else:
        header = next(
            csv.reader(
                [first_line],
                delimiter=delimiter
            )
        )

    second_line = fh.readline()

    if not second_line:
        raise RuntimeError(
            "UMI matrix contains header only"
        )

    if delimiter == " ":
        first_data = second_line.strip().split()
    else:
        first_data = next(
            csv.reader(
                [second_line],
                delimiter=delimiter
            )
        )


print(
    "Header columns:",
    len(header)
)

print(
    "First data row columns:",
    len(first_data)
)


# ------------------------------------------------------------
# Figure out whether header contains a gene-name column
#
# GSE160269 format (Case B):
# barcode1 barcode2 barcode3 ...
# GENE1 0 1 0 ...
#
# header = barcodes
# first_data = [gene_name, val1, val2, ...]
# ------------------------------------------------------------

if len(first_data) == len(header) + 1:

    barcodes = header
    header_has_gene_col = False

    print(
        "Format: barcodes in header, gene name in data rows"
    )

elif len(first_data) == len(header):

    barcodes = header[1:]
    header_has_gene_col = True

    print(
        "Format: gene name column in header"
    )

else:

    raise RuntimeError(
        "Could not determine matrix orientation/header format.\n"
        f"Header={len(header)}, first data row={len(first_data)}"
    )


n_cells = len(barcodes)

print(
    "Cells detected:",
    f"{n_cells:,}"
)

if n_cells == 0:
    raise RuntimeError(
        "No cell barcodes detected."
    )


# ------------------------------------------------------------
# Save barcodes
# ------------------------------------------------------------

with gzip.open(
    barcode_out,
    "wt",
    encoding="utf-8"
) as f:

    for bc in barcodes:
        f.write(
            str(bc).strip() + "\n"
        )


# ------------------------------------------------------------
# Stream through gene rows and write ONLY non-zero values.
# Peak memory stays small.
# ------------------------------------------------------------

genes = []
nnz = 0
gene_index = 0

start = time.time()


def process_row(row, entries):

    global gene_index, nnz

    if not row:
        return

    expected = n_cells + 1

    if len(row) != expected:

        raise RuntimeError(
            f"Row {gene_index + 1}: "
            f"expected {expected} columns, "
            f"found {len(row)}"
        )

    gene = row[0].strip()

    gene_index += 1

    if not gene:
        gene = f"Gene_{gene_index}"

    genes.append(gene)

    values = row[1:]

    for j, raw in enumerate(
        values,
        start=1
    ):

        raw = raw.strip()

        if (
            raw == "" or
            raw == "0" or
            raw == "0.0"
        ):
            continue

        try:
            value = float(raw)
        except ValueError:
            raise RuntimeError(
                f"Non-numeric value at "
                f"gene {gene}, cell {j}: {raw}"
            )

        if value == 0:
            continue

        # MatrixMarket: row col value
        entries.write(
            f"{gene_index} {j} {raw}\n"
        )

        nnz += 1


with open(
    tmp_entries,
    "w",
    encoding="utf-8"
) as entries:

    # process first data row
    process_row(
        first_data,
        entries
    )

    with open_text(source) as fh:

        # Skip header
        next(fh)

        # Skip first data row because already processed
        next(fh)

        if delimiter == " ":
            # Space-delimited: split each line on whitespace
            for line in fh:
                row = line.strip().split()
                if row:
                    process_row(
                        row,
                        entries
                    )

                    if gene_index % 500 == 0:
                        elapsed = (
                            time.time() - start
                        )
                        print(
                            f"[{label}] "
                            f"genes={gene_index:,} "
                            f"nnz={nnz:,} "
                            f"elapsed={elapsed/60:.1f} min",
                            flush=True
                        )
        else:
            reader = csv.reader(
                fh,
                delimiter=delimiter
            )

            for row in reader:

                process_row(
                    row,
                    entries
                )

                if gene_index % 500 == 0:

                    elapsed = (
                        time.time() - start
                    )

                    print(
                        f"[{label}] "
                        f"genes={gene_index:,} "
                        f"nnz={nnz:,} "
                        f"elapsed={elapsed/60:.1f} min",
                        flush=True
                    )


n_genes = gene_index

print()
print(
    "Genes:",
    f"{n_genes:,}"
)

print(
    "Cells:",
    f"{n_cells:,}"
)

print(
    "Non-zero entries:",
    f"{nnz:,}"
)


# ------------------------------------------------------------
# Save features
# ------------------------------------------------------------

with gzip.open(
    feature_out,
    "wt",
    encoding="utf-8"
) as f:

    for gene in genes:
        f.write(
            gene + "\n"
        )


# ------------------------------------------------------------
# Build final gzipped MatrixMarket
# ------------------------------------------------------------

print(
    "Writing MatrixMarket..."
)

with gzip.open(
    matrix_out,
    "wt",
    encoding="utf-8"
) as out:

    out.write(
        "%%MatrixMarket matrix coordinate real general\n"
    )

    out.write(
        "% converted from GSE160269 wide UMI matrix\n"
    )

    out.write(
        f"{n_genes} {n_cells} {nnz}\n"
    )

    with open(
        tmp_entries,
        "r",
        encoding="utf-8"
    ) as entries:

        shutil.copyfileobj(
            entries,
            out
        )


os.remove(
    tmp_entries
)

with open(
    os.path.join(
        outdir,
        ".conversion_complete"
    ),
    "w"
) as f:

    f.write(
        f"genes={n_genes}\n"
    )

    f.write(
        f"cells={n_cells}\n"
    )

    f.write(
        f"nnz={nnz}\n"
    )


print(
    "✓ Conversion finished:",
    label
)
PY

done


echo ""
echo "============================================================"
echo "3. BUILD ONE SEURAT OBJECT PER MATRIX"
echo "============================================================"
echo ""


# ============================================================
# 3. Build objects sequentially
# ============================================================

cat > "$PROJECT/07_scripts/12B_GSE160269_build_parts.R" <<'RSCRIPT'

options(
  stringsAsFactors = FALSE
)

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
})

proc_dir <- "02_processed/GSE160269/10x_parts"
obj_dir  <- "03_objects/GSE160269/parts"

dir.create(
  obj_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

manifest <- fread(
  "00_metadata/06_GSE160269_UMI_matrix_manifest.tsv"
)

audit <- list()


for (i in seq_len(nrow(manifest))) {

  label <- manifest$label[i]

  cat(
    "\n============================================================\n"
  )

  cat(
    "Building:",
    label,
    "\n"
  )

  cat(
    "============================================================\n"
  )

  part_dir <- file.path(
    proc_dir,
    label
  )

  out_file <- file.path(
    obj_dir,
    paste0(
      "GSE160269_",
      label,
      ".rds"
    )
  )


  # ----------------------------------------------------------
  # Resume support
  # ----------------------------------------------------------

  if (file.exists(out_file)) {

    cat(
      "✓ Object already exists — skipping.\n"
    )

    info <- file.info(
      out_file
    )

    audit[[label]] <- data.table(
      label = label,
      object_file = out_file,
      status = "existing",
      object_size_MB =
        round(
          info$size / 1024^2,
          2
        )
    )

    next
  }


  # ----------------------------------------------------------
  # Read sparse matrix
  # ----------------------------------------------------------

  counts <- ReadMtx(
    mtx = file.path(
      part_dir,
      "matrix.mtx.gz"
    ),

    cells = file.path(
      part_dir,
      "barcodes.tsv.gz"
    ),

    features = file.path(
      part_dir,
      "features.tsv.gz"
    ),

    cell.column = 1,
    feature.column = 1,
    unique.features = TRUE
  )


  original_barcode <- colnames(
    counts
  )


  # ----------------------------------------------------------
  # IMPORTANT:
  # Make globally unique cell names BEFORE any future merge.
  #
  # This prevents the exact problem that broke the old script.
  # ----------------------------------------------------------

  globally_unique_name <- paste(
    label,
    original_barcode,
    sep = "__"
  )

  colnames(
    counts
  ) <- globally_unique_name


  # ----------------------------------------------------------
  # Build Seurat object
  # ----------------------------------------------------------

  obj <- CreateSeuratObject(
    counts = counts,
    project = "GSE160269",
    min.cells = 0,
    min.features = 0
  )


  obj$dataset <- "GSE160269"

  obj$source_matrix <- label

  obj$source_celltype <- label

  obj$original_barcode <-
    original_barcode


  cat(
    "Genes:",
    nrow(obj),
    "\n"
  )

  cat(
    "Cells:",
    ncol(obj),
    "\n"
  )


  # ----------------------------------------------------------
  # Save immediately
  # ----------------------------------------------------------

  saveRDS(
    obj,
    out_file,
    compress = FALSE
  )


  info <- file.info(
    out_file
  )


  audit[[label]] <- data.table(
    label = label,
    genes = nrow(obj),
    cells = ncol(obj),
    object_file = out_file,
    status = "built",
    object_size_MB =
      round(
        info$size / 1024^2,
        2
      )
  )


  cat(
    "✓ Saved:",
    out_file,
    "\n"
  )


  # ----------------------------------------------------------
  # Explicitly release RAM before next matrix
  # ----------------------------------------------------------

  rm(
    obj,
    counts
  )

  invisible(
    gc()
  )
}


# ============================================================
# Audit
# ============================================================

audit_dt <- rbindlist(
  audit,
  fill = TRUE
)

fwrite(
  audit_dt,
  "00_metadata/06_GSE160269_part_object_audit.csv",
  bom = TRUE
)

cat(
  "\n============================================================\n"
)

cat(
  "GSE160269 PART OBJECT AUDIT\n"
)

cat(
  "============================================================\n\n"
)

print(
  audit_dt
)


cat(
  "\nObjects present:\n"
)

objects <- list.files(
  obj_dir,
  pattern = "\\.rds$",
  full.names = TRUE
)

for (f in objects) {

  info <- file.info(f)

  cat(
    sprintf(
      "  %-50s %8.2f MB\n",
      basename(f),
      info$size / 1024^2
    )
  )
}


cat(
  "\n============================================================\n"
)

cat(
  "STEP 13B FINISHED ✓\n"
)

cat(
  "============================================================\n\n"
)

cat(
  "IMPORTANT:\n"
)

cat(
  "  GSE160269 parts were built separately.\n"
)

cat(
  "  They were NOT merged.\n"
)

cat(
  "  Metadata matching was NOT attempted yet.\n"
)

cat(
  "  QC was NOT run.\n"
)

cat(
  "  Normalization was NOT run.\n"
)

cat(
  "  Integration was NOT run.\n"
)

cat(
  "  STEP 14 was NOT run.\n\n"
)

cat(
  "STOP HERE and review the audit before continuing.\n"
)

RSCRIPT


Rscript \
  07_scripts/12B_GSE160269_build_parts.R \
  2>&1 |
tee \
  08_logs/12B_GSE160269_build_parts.log


echo ""
echo "============================================================"
echo "FINAL FILE CHECK"
echo "============================================================"
echo ""

find \
  "$OBJ_DIR" \
  -maxdepth 1 \
  -type f \
  -name "*.rds" \
  -exec ls -lh {} \;

echo ""
echo "============================================================"
echo "STEP 13B COMPLETE"
echo ""
echo "STOP."
echo "DO NOT RUN STEP 14."
echo "============================================================"
