#!/bin/bash

set -euo pipefail

PROJECT="$HOME/Downloads/ESCC Neoadjuvant"

META_DIR="$PROJECT/00_metadata/OMIX005710"
RAW_DIR="$PROJECT/01_raw/OMIX005710"
PAIRED_DIR="$RAW_DIR/paired_tumor"
LOG_DIR="$PROJECT/08_logs"

PLAN="$META_DIR/STEP22A_TARGETED_TUMOR_DOWNLOAD_PLAN.csv"

cd "$PROJECT"

mkdir -p "$PAIRED_DIR"
mkdir -p "$LOG_DIR"

echo ""
echo "========================================================================"
echo "STEP22A-3A: DOWNLOAD PRIMARY PAIRED TUMOR FILES"
echo "========================================================================"
echo ""

# ======================================================================
# 1. Guards
# ======================================================================

echo "===== 1. PROJECT / INPUT CHECK ====="
echo ""

if [ ! -f "$PLAN" ]; then
    echo "ERROR: missing download plan:"
    echo "$PLAN"
    exit 1
fi

echo "✓ Download plan found:"
echo "$PLAN"
echo ""

# ======================================================================
# 2. Build paired-only target list automatically
# ======================================================================

echo "===== 2. BUILD PAIRED-ONLY DOWNLOAD TARGETS ====="
echo ""

python3 - <<'PY'
import csv
import os
import sys

plan = (
    "00_metadata/OMIX005710/"
    "STEP22A_TARGETED_TUMOR_DOWNLOAD_PLAN.csv"
)

out = (
    "00_metadata/OMIX005710/"
    "STEP22A_03A_PAIRED_DOWNLOAD_TARGETS.tsv"
)

expected_patients = {
    "P6",
    "P9",
    "P15",
    "P16",
    "P18",
    "P19",
}

with open(
    plan,
    encoding="utf-8-sig",
    newline=""
) as f:
    rows = list(csv.DictReader(f))

selected = []

for r in rows:

    patient = (r.get("patient_id") or "").strip()
    timepoint = (r.get("timepoint") or "").strip()
    sample = (r.get("sample_name") or "").strip()
    hrs = (r.get("hrs_accession") or "").strip()
    file_id = (r.get("omix_file_id") or "").strip()
    size_mb = (r.get("listed_size_MB") or "").strip()
    url = (r.get("download_url_https") or "").strip()
    priority = (r.get("priority") or "").strip()

    if patient not in expected_patients:
        continue

    if priority != "1_PRIMARY_PAIRED":
        raise SystemExit(
            f"ERROR: {sample} belongs to paired patient "
            f"but priority is {priority}"
        )

    if timepoint not in {"Before", "After"}:
        raise SystemExit(
            f"ERROR: invalid timepoint for {sample}: "
            f"{timepoint}"
        )

    if not all([
        patient,
        sample,
        hrs,
        file_id,
        url
    ]):
        raise SystemExit(
            f"ERROR: incomplete target row for {sample}"
        )

    if not url.startswith("https://download.cncb.ac.cn/"):
        raise SystemExit(
            "ERROR: unexpected download host:\n"
            + url
        )

    filename = file_id + ".tar"

    selected.append({
        "patient_id": patient,
        "timepoint": timepoint,
        "sample_name": sample,
        "hrs_accession": hrs,
        "file_id": file_id,
        "filename": filename,
        "size_mb": size_mb,
        "url": url,
    })


def pkey(x):
    return int(
        "".join(
            c for c in x["patient_id"]
            if c.isdigit()
        )
    )


selected.sort(
    key=lambda x: (
        pkey(x),
        0 if x["timepoint"] == "Before" else 1
    )
)

if len(selected) != 12:
    raise SystemExit(
        f"ERROR: expected 12 paired tumor files, "
        f"found {len(selected)}"
    )

patients = {
    x["patient_id"]
    for x in selected
}

if patients != expected_patients:
    raise SystemExit(
        "ERROR: paired patient set mismatch:\n"
        + str(sorted(patients))
    )

for p in sorted(
    expected_patients,
    key=lambda x: int(x[1:])
):

    tp = {
        x["timepoint"]
        for x in selected
        if x["patient_id"] == p
    }

    if tp != {"Before", "After"}:
        raise SystemExit(
            f"ERROR: {p} does not have both "
            "Before and After."
        )

with open(
    out,
    "w",
    encoding="utf-8",
    newline=""
) as f:

    fields = [
        "patient_id",
        "timepoint",
        "sample_name",
        "hrs_accession",
        "file_id",
        "filename",
        "size_mb",
        "url",
    ]

    writer = csv.DictWriter(
        f,
        fieldnames=fields,
        delimiter="\t"
    )

    writer.writeheader()
    writer.writerows(selected)


total = 0.0

for x in selected:

    try:
        total += float(x["size_mb"])
    except Exception:
        pass


print("Paired patients:")
print(
    ", ".join(
        sorted(
            patients,
            key=lambda x: int(x[1:])
        )
    )
)

print()
print("Download files:", len(selected))
print(
    "Listed total size:",
    f"{total:.2f} MB "
    f"({total / 1024:.3f} GiB)"
)

print()
print("Targets:")

for x in selected:
    print(
        f"  {x['patient_id']:>3s}  "
        f"{x['timepoint']:<6s}  "
        f"{x['sample_name']:<12s}  "
        f"{x['filename']}"
    )

print()
print("Saved:")
print(out)
PY

echo ""

# ======================================================================
# 3. Disk-space guard
# ======================================================================

echo "===== 3. DISK SPACE CHECK ====="
echo ""

df -h "$PROJECT"
echo ""

FREE_KB=$(df -Pk "$PROJECT" | awk 'NR==2 {print $4}')

# Keep at least 20 GiB free before starting.
MIN_FREE_KB=$((20 * 1024 * 1024))

if [ "$FREE_KB" -lt "$MIN_FREE_KB" ]; then
    echo "ERROR: less than 20 GiB free."
    echo "STEP22A_INSUFFICIENT_DISK_SPACE"
    exit 1
fi

echo "✓ More than 20 GiB free."
echo ""

# ======================================================================
# 4. Download with resume
# ======================================================================

echo "===== 4. DOWNLOAD PAIRED TUMOR FILES ====="
echo ""

TARGETS="$META_DIR/STEP22A_03A_PAIRED_DOWNLOAD_TARGETS.tsv"

tail -n +2 "$TARGETS" | tr -d '\r' |
while IFS=$'\t' read -r \
    patient \
    timepoint \
    sample \
    hrs \
    file_id \
    filename \
    size_mb \
    url

do

    final="$PAIRED_DIR/$filename"
    partial="$final.part"

    echo ""
    echo "------------------------------------------------------------------------"
    echo "PATIENT   : $patient"
    echo "TIMEPOINT : $timepoint"
    echo "SAMPLE    : $sample"
    echo "HRS       : $hrs"
    echo "FILE      : $filename"
    echo "EXPECTED  : ${size_mb} MB"
    echo "------------------------------------------------------------------------"
    echo ""

    if [ -f "$final" ]; then

        echo "[SKIP] Final file already exists."
        echo "$final"

    else

        echo "[DOWNLOAD / RESUME]"

        curl \
          -L \
          --fail \
          --show-error \
          --retry 8 \
          --retry-delay 10 \
          --retry-all-errors \
          --connect-timeout 30 \
          -C - \
          -o "$partial" \
          "$url"

        mv "$partial" "$final"

        echo ""
        echo "[DONE]"
        echo "$final"

    fi

done

echo ""

# ======================================================================
# 5. Check no .part remains
# ======================================================================

echo "===== 5. PARTIAL-FILE CHECK ====="
echo ""

PART_COUNT=$(find "$PAIRED_DIR" \
    -maxdepth 1 \
    -type f \
    -name '*.part' \
    | wc -l \
    | tr -d ' ')

echo "Remaining .part files: $PART_COUNT"

if [ "$PART_COUNT" -ne 0 ]; then
    echo ""
    echo "DOWNLOAD_NOT_COMPLETE"
    echo "Do NOT delete .part files."
    exit 1
fi

echo "✓ No .part files remain."
echo ""

# ======================================================================
# 6. Count downloaded TAR files
# ======================================================================

echo "===== 6. DOWNLOAD COUNT ====="
echo ""

TAR_COUNT=$(find "$PAIRED_DIR" \
    -maxdepth 1 \
    -type f \
    -name '*.tar' \
    | wc -l \
    | tr -d ' ')

echo "Downloaded TAR files: $TAR_COUNT"

if [ "$TAR_COUNT" -ne 12 ]; then
    echo ""
    echo "ERROR: expected exactly 12 TAR files."
    exit 1
fi

echo "✓ Exactly 12 paired tumor TAR files present."
echo ""

# ======================================================================
# 7. TAR readability check
# ======================================================================

echo "===== 7. TAR READABILITY CHECK ====="
echo ""

FAIL=0

for f in "$PAIRED_DIR"/*.tar
do

    name=$(basename "$f")

    printf "%-35s " "$name"

    if tar -tf "$f" >/dev/null 2>&1; then
        echo "✓ READABLE"
    else
        echo "✗ FAILED"
        FAIL=$((FAIL + 1))
    fi

done

if [ "$FAIL" -ne 0 ]; then
    echo ""
    echo "ERROR: unreadable TAR archive(s): $FAIL"
    exit 1
fi

echo ""
echo "✓ All 12 TAR archives are readable."
echo ""

# ======================================================================
# 8. Local size + hashes + provenance
# ======================================================================

echo "===== 8. FILE SIZE + SHA256 + MD5 ====="
echo ""

python3 - <<'PY'
import csv
import hashlib
import os
from datetime import datetime

project = os.getcwd()

targets_file = (
    "00_metadata/OMIX005710/"
    "STEP22A_03A_PAIRED_DOWNLOAD_TARGETS.tsv"
)

data_dir = (
    "01_raw/OMIX005710/"
    "paired_tumor"
)

outfile = (
    "00_metadata/OMIX005710/"
    "STEP22A_03A_DOWNLOAD_PROVENANCE.csv"
)

with open(
    targets_file,
    encoding="utf-8",
    newline=""
) as f:

    targets = list(
        csv.DictReader(
            f,
            delimiter="\t"
        )
    )

rows = []

for t in targets:

    filename = t["filename"]

    path = os.path.join(
        data_dir,
        filename
    )

    if not os.path.exists(path):
        raise SystemExit(
            "ERROR missing file: " + path
        )

    size_bytes = os.path.getsize(path)
    actual_mb = size_bytes / 1024**2

    sha = hashlib.sha256()
    md5 = hashlib.md5()

    print(
        "Hashing:",
        t["patient_id"],
        t["timepoint"],
        filename
    )

    with open(path, "rb") as fh:

        while True:

            chunk = fh.read(
                8 * 1024 * 1024
            )

            if not chunk:
                break

            sha.update(chunk)
            md5.update(chunk)

    expected_mb = None

    try:
        expected_mb = float(
            t["size_mb"]
        )
    except Exception:
        pass

    if expected_mb is None:

        size_status = (
            "EXPECTED_SIZE_UNAVAILABLE"
        )

        difference_mb = ""

    else:

        difference = (
            actual_mb - expected_mb
        )

        difference_mb = round(
            difference,
            3
        )

        # OMIX displayed sizes are rounded.
        tolerance = max(
            1.0,
            expected_mb * 0.02
        )

        size_status = (
            "PASS"
            if abs(difference) <= tolerance
            else "REVIEW"
        )

    rows.append({

        "patient_id":
            t["patient_id"],

        "timepoint":
            t["timepoint"],

        "sample_name":
            t["sample_name"],

        "hrs_accession":
            t["hrs_accession"],

        "omix_file_id":
            t["file_id"],

        "filename":
            filename,

        "url":
            t["url"],

        "expected_size_MB":
            t["size_mb"],

        "actual_size_bytes":
            size_bytes,

        "actual_size_MB":
            round(actual_mb, 3),

        "size_difference_MB":
            difference_mb,

        "size_check":
            size_status,

        "md5":
            md5.hexdigest(),

        "sha256":
            sha.hexdigest(),

        "download_checked_at":
            datetime.now().isoformat(
                timespec="seconds"
            ),

        "tar_readable":
            "YES",

        "status":
            "DOWNLOADED_VERIFIED"
    })


with open(
    outfile,
    "w",
    newline="",
    encoding="utf-8-sig"
) as f:

    writer = csv.DictWriter(
        f,
        fieldnames=list(
            rows[0].keys()
        )
    )

    writer.writeheader()
    writer.writerows(rows)


print()
print("Verified files:", len(rows))

print(
    "Total actual size:",
    round(
        sum(
            x["actual_size_bytes"]
            for x in rows
        ) / 1024**3,
        3
    ),
    "GiB"
)

review = [
    x for x in rows
    if x["size_check"] == "REVIEW"
]

print(
    "Size checks requiring review:",
    len(review)
)

print()
print("Saved:")
print(outfile)
PY

echo ""

# ======================================================================
# 9. Inspect TAR structure WITHOUT extracting
# ======================================================================

echo "===== 9. TAR CONTENT AUDIT ====="
echo ""

python3 - <<'PY'
import csv
import os
import tarfile
from collections import Counter

targets_file = (
    "00_metadata/OMIX005710/"
    "STEP22A_03A_PAIRED_DOWNLOAD_TARGETS.tsv"
)

data_dir = (
    "01_raw/OMIX005710/"
    "paired_tumor"
)

outfile = (
    "00_metadata/OMIX005710/"
    "STEP22A_03A_TAR_CONTENT_AUDIT.csv"
)

with open(
    targets_file,
    encoding="utf-8",
    newline=""
) as f:

    targets = list(
        csv.DictReader(
            f,
            delimiter="\t"
        )
    )

audit = []

for t in targets:

    path = os.path.join(
        data_dir,
        t["filename"]
    )

    with tarfile.open(
        path,
        "r:*"
    ) as tar:

        members = [
            m
            for m in tar.getmembers()
            if m.isfile()
        ]

        names = [
            m.name
            for m in members
        ]

    lower = [
        x.lower()
        for x in names
    ]

    matrix = [
        x for x in names
        if (
            x.lower().endswith(
                ".mtx"
            )
            or
            x.lower().endswith(
                ".mtx.gz"
            )
        )
    ]

    barcodes = [
        x for x in names
        if "barcode" in x.lower()
    ]

    features = [
        x for x in names
        if (
            "feature" in x.lower()
            or
            "gene" in x.lower()
        )
        and (
            x.lower().endswith(
                ".tsv"
            )
            or
            x.lower().endswith(
                ".tsv.gz"
            )
            or
            x.lower().endswith(
                ".txt"
            )
            or
            x.lower().endswith(
                ".txt.gz"
            )
        )
    ]

    annotation = [
        x for x in names
        if any(
            k in x.lower()
            for k in [
                "annot",
                "celltype",
                "cell_type",
                "metadata",
                "cluster"
            ]
        )
    ]

    ext_counter = Counter()

    for name in lower:

        if name.endswith(".mtx.gz"):
            ext = ".mtx.gz"

        elif name.endswith(".tsv.gz"):
            ext = ".tsv.gz"

        elif name.endswith(".csv.gz"):
            ext = ".csv.gz"

        elif name.endswith(".txt.gz"):
            ext = ".txt.gz"

        elif name.endswith(".mtx"):
            ext = ".mtx"

        elif name.endswith(".tsv"):
            ext = ".tsv"

        elif name.endswith(".csv"):
            ext = ".csv"

        elif name.endswith(".txt"):
            ext = ".txt"

        else:
            ext = (
                os.path.splitext(name)[1]
                or "other"
            )

        ext_counter[ext] += 1

    if (
        matrix
        and
        barcodes
        and
        features
    ):
        format_status = (
            "10X_STYLE_COMPONENTS_DETECTED"
        )

    else:
        format_status = (
            "NONSTANDARD_OR_REVIEW_REQUIRED"
        )

    audit.append({

        "patient_id":
            t["patient_id"],

        "timepoint":
            t["timepoint"],

        "sample_name":
            t["sample_name"],

        "hrs_accession":
            t["hrs_accession"],

        "filename":
            t["filename"],

        "internal_file_count":
            len(names),

        "matrix_files":
            " | ".join(matrix),

        "barcode_files":
            " | ".join(barcodes),

        "feature_or_gene_files":
            " | ".join(features),

        "annotation_files":
            " | ".join(annotation),

        "extension_summary":
            "; ".join(
                f"{k}:{v}"
                for k, v
                in sorted(
                    ext_counter.items()
                )
            ),

        "format_status":
            format_status
    })


with open(
    outfile,
    "w",
    newline="",
    encoding="utf-8-sig"
) as f:

    writer = csv.DictWriter(
        f,
        fieldnames=list(
            audit[0].keys()
        )
    )

    writer.writeheader()
    writer.writerows(audit)


print(
    f"{'PATIENT':<8}"
    f"{'TIME':<8}"
    f"{'SAMPLE':<12}"
    f"{'FILES':<8}"
    f"FORMAT"
)

print("-" * 75)

for x in audit:

    print(
        f"{x['patient_id']:<8}"
        f"{x['timepoint']:<8}"
        f"{x['sample_name']:<12}"
        f"{x['internal_file_count']:<8}"
        f"{x['format_status']}"
    )

print()

print(
    "10x-style archives:",
    sum(
        x["format_status"]
        ==
        "10X_STYLE_COMPONENTS_DETECTED"
        for x in audit
    ),
    "/",
    len(audit)
)

print()

print("Annotation-like files detected:")

annotation_hits = 0

for x in audit:

    if x["annotation_files"]:

        annotation_hits += 1

        print(
            " ",
            x["sample_name"],
            "->",
            x["annotation_files"]
        )

if annotation_hits == 0:
    print(
        "  <none inside paired TAR files>"
    )

print()
print("Saved:")
print(outfile)
PY

echo ""

# ======================================================================
# 10. Show first TAR file contents as example
# ======================================================================

echo "===== 10. EXAMPLE TAR CONTENT ====="
echo ""

FIRST_TAR=$(find "$PAIRED_DIR" \
    -maxdepth 1 \
    -type f \
    -name '*.tar' \
    | sort \
    | head -1)

echo "Example archive:"
echo "$FIRST_TAR"
echo ""

tar -tf "$FIRST_TAR" | head -60

echo ""

# ======================================================================
# 11. Final disk space
# ======================================================================

echo "===== 11. DISK SPACE AFTER DOWNLOAD ====="
echo ""

df -h "$PROJECT"

echo ""

# ======================================================================
# FINAL
# ======================================================================

echo "========================================================================"
echo "STEP22A-3A COMPLETE"
echo "========================================================================"
echo ""

echo "PAIRED PATIENTS              : 6"
echo "PAIRED TUMOR FILES           : 12"
echo "DOWNLOAD MODE                : PRIMARY PAIRED ONLY"
echo ""
echo "UNPAIRED TUMOR FILES         : NOT DOWNLOADED"
echo "ADJACENT NORMAL FILES        : NOT DOWNLOADED"
echo ""
echo "TAR EXTRACTION               : NOT PERFORMED"
echo "TAR READABILITY              : VERIFIED"
echo "LOCAL MD5 / SHA256           : RECORDED"
echo "TAR INTERNAL STRUCTURE       : AUDITED"
echo ""
echo "FROZEN STEP19-21 MODIFIED    : NO"
echo "MANUSCRIPT MODIFIED          : NO"
echo ""
echo "READY FOR STEP22A-3B:"
echo "PAIRED TAR EXTRACTION + MATRIX/ANNOTATION AUDIT"
echo ""
echo "STOP HERE."
echo "DO NOT MANUALLY EXTRACT TAR FILES."
echo ""
