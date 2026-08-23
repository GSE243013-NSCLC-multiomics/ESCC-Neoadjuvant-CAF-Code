#!/bin/bash

set -e

PROJECT="$HOME/Downloads/ESCC Neoadjuvant"
RAW_DIR="$PROJECT/01_raw/GEO/GSE160269"
OUT_DIR="$PROJECT/06_tables/GSE160269/metadata"

mkdir -p "$OUT_DIR"

cd "$PROJECT"

echo "============================================================"
echo "STEP 13C2: FIND REAL GSE160269 CELL METADATA"
echo ""
echo "NO MERGE"
echo "NO OBJECT MODIFICATION"
echo "NO QC"
echo "NO NORMALIZATION"
echo "NO INTEGRATION"
echo "============================================================"

echo ""
echo "===== 1. ALL 13 RAW FILES ====="
echo ""

find "$RAW_DIR" \
  -maxdepth 2 \
  -type f \
  -exec ls -lh {} \; | sort


echo ""
echo "===== 2. CANDIDATE CELL METADATA FILES ====="
echo ""

python3 - <<'PY'
import os
import re

root = "01_raw/GEO/GSE160269"

all_files = []

for dp, _, fs in os.walk(root):
    for fn in fs:
        all_files.append(os.path.join(dp, fn))

candidates = []

for path in sorted(all_files):

    fn = os.path.basename(path)
    low = fn.lower()

    # text-like files only
    if not re.search(r"\.(csv|tsv|txt)(\.gz)?$", low):
        continue

    # exclude matrices and special datasets
    if "umi" in low:
        continue

    if "tcr" in low:
        continue

    if "cd45" in low and "meta" not in low and "annot" not in low:
        continue

    # probable metadata / annotation
    if any(
        x in low
        for x in [
            "meta",
            "annot",
            "celltype",
            "cell_type",
            "cluster"
        ]
    ):
        candidates.append(path)

print(f"Detected {len(candidates)} probable non-TCR metadata files:\n")

for i, p in enumerate(candidates, 1):

    size = os.path.getsize(p) / 1024**2

    print(
        f"{i:02d}. {os.path.basename(p)}"
        f"  ({size:.2f} MB)"
    )

with open(
    "06_tables/GSE160269/metadata/"
    "GSE160269_real_metadata_candidates.txt",
    "w",
    encoding="utf-8"
) as f:

    for p in candidates:
        f.write(p + "\n")
PY


echo ""
echo "===== 3. INSPECT HEADERS AND FIRST 5 ROWS ====="
echo ""

python3 - <<'PY'
import os
import gzip
import csv
import re

candidate_file = (
    "06_tables/GSE160269/metadata/"
    "GSE160269_real_metadata_candidates.txt"
)

if not os.path.exists(candidate_file):
    raise SystemExit("Candidate list missing")

with open(candidate_file, encoding="utf-8") as f:
    files = [
        x.strip()
        for x in f
        if x.strip()
    ]


def open_text(path):

    if path.lower().endswith(".gz"):
        return gzip.open(
            path,
            "rt",
            encoding="utf-8-sig",
            errors="replace",
            newline=""
        )

    return open(
        path,
        "r",
        encoding="utf-8-sig",
        errors="replace",
        newline=""
    )


for path in files:

    print()
    print("=" * 90)
    print(os.path.basename(path))
    print("=" * 90)

    with open_text(path) as fh:

        first = fh.readline()

        if not first:
            print("EMPTY FILE")
            continue

        # delimiter guess
        if first.count("\t") >= first.count(","):
            delim = "\t"
        else:
            delim = ","

        header = next(
            csv.reader(
                [first],
                delimiter=delim
            )
        )

        print(
            f"Delimiter: {repr(delim)}"
        )

        print(
            f"Columns: {len(header)}"
        )

        print()
        print("COLUMN NAMES:")

        for i, col in enumerate(header, 1):
            print(
                f"  {i:02d}. {col}"
            )

        print()
        print("FIRST 5 ROWS:")

        reader = csv.reader(
            fh,
            delimiter=delim
        )

        for i, row in enumerate(reader, 1):

            if i > 5:
                break

            # avoid giant lines
            preview = row[:20]

            print(
                f"\nROW {i}:"
            )

            print(preview)
PY


echo ""
echo "===== 4. CHECK THE 8 SEURAT PART CELL-ID FORMAT ====="
echo ""

Rscript - <<'RSCRIPT'

suppressPackageStartupMessages({
  library(Seurat)
})

obj_dir <- "03_objects/GSE160269/parts"

files <- sort(
  list.files(
    obj_dir,
    pattern = "\\.rds$",
    full.names = TRUE
  )
)

for (f in files) {

  obj <- readRDS(f)

  cat(
    "\n============================================================\n"
  )

  cat(
    basename(f),
    "\n"
  )

  cat(
    "Cells:",
    ncol(obj),
    "\n"
  )

  cat(
    "First 5 current cell IDs:\n"
  )

  print(
    head(
      colnames(obj),
      5
    )
  )

  cat(
    "First 5 original barcodes:\n"
  )

  print(
    head(
      obj$original_barcode,
      5
    )
  )

  rm(obj)

  gc(verbose = FALSE)
}

RSCRIPT


echo ""
echo "============================================================"
echo "STEP 13C2 FINISHED"
echo ""
echo "STOP HERE."
echo "The previous 32.9% result was TCR-specific."
echo "No objects were modified."
echo "STEP 14 was NOT run."
echo "============================================================"
