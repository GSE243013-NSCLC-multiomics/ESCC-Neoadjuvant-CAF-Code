#!/bin/bash
set -e

PROJECT="$HOME/Downloads/ESCC Neoadjuvant"
DATA_DIR="$PROJECT/01_raw/GEO/GSE221561"
META_DIR="$PROJECT/00_metadata"
LOG_DIR="$PROJECT/08_logs"

cd "$PROJECT"
mkdir -p "$DATA_DIR" "$META_DIR" "$LOG_DIR"

echo ""
echo "============================================================"
echo "STEP 7: DOWNLOAD GSE221561"
echo "============================================================"
echo ""

# 1. Check disk space
echo "===== 1. Disk space ====="
df -h "$PROJECT"
echo ""

# 2. Extract download URLs from plan
echo "===== 2. Preparing download targets ====="
python3 - <<'PY'
import csv, os, sys

plan = "00_metadata/01_GEO_download_plan.csv"
out = "00_metadata/03_GSE221561_download_targets.tsv"

if not os.path.exists(plan):
    sys.exit(f"ERROR: cannot find {plan}")

targets = []
with open(plan, encoding="utf-8-sig") as f:
    reader = csv.DictReader(f)
    for row in reader:
        if row.get("GSE") == "GSE221561":
            filename = (row.get("filename") or "").strip()
            url = (row.get("source_url") or "").strip()
            if filename and url:
                targets.append((filename, url))

if len(targets) != 2:
    print(f"WARNING: expected 2 files, found {len(targets)}")

with open(out, "w", encoding="utf-8") as f:
    for filename, url in targets:
        f.write(f"{filename}\t{url}\n")

print("Download targets:")
for filename, url in targets:
    print(f"  {filename}")
    print(f"  {url}")
    print()
print("Saved:", out)
PY

echo ""

# 3. Download
echo "===== 3. Downloading GSE221561 ====="
echo "Large file is about 916 MB."
echo ""

while IFS=$'\t' read -r filename url; do
    destination="$DATA_DIR/$filename"
    partial="$destination.part"

    echo "------------------------------------------------------------"
    echo "FILE: $filename"
    echo "URL : $url"
    echo "------------------------------------------------------------"

    if [ -f "$destination" ]; then
        echo "[SKIP] File already exists: $destination"
    else
        if [ -f "$partial" ]; then
            echo "[RESUME] Continuing interrupted download..."
            curl -L --fail --retry 5 --retry-delay 5 --connect-timeout 30 -C - -o "$partial" "$url"
        else
            echo "[DOWNLOAD] Starting..."
            curl -L --fail --retry 5 --retry-delay 5 --connect-timeout 30 -o "$partial" "$url"
        fi
        mv "$partial" "$destination"
        echo "[DONE] $filename"
    fi
    echo ""
done < "$META_DIR/03_GSE221561_download_targets.tsv"

# 4. File verification
echo ""
echo "===== 4. Downloaded files ====="
ls -lh "$DATA_DIR"
echo ""
file "$DATA_DIR"/*
echo ""

# 5. MD5 manifest
echo "===== 5. Calculating MD5 ====="
python3 - <<'PY'
import hashlib, csv, os
from datetime import datetime

data_dir = "01_raw/GEO/GSE221561"
outfile = "00_metadata/03_GSE221561_download_manifest.csv"

files = sorted([f for f in os.listdir(data_dir) if os.path.isfile(os.path.join(data_dir, f)) and not f.endswith(".part")])
rows = []

for filename in files:
    path = os.path.join(data_dir, filename)
    print("Calculating MD5:", filename)
    md5 = hashlib.md5()
    with open(path, "rb") as fh:
        while True:
            chunk = fh.read(8 * 1024 * 1024)
            if not chunk: break
            md5.update(chunk)
    size = os.path.getsize(path)
    rows.append({"GSE": "GSE221561", "filename": filename, "size_bytes": size, "size_MB": round(size / 1024**2, 2), "md5": md5.hexdigest(), "download_verified_at": datetime.now().isoformat(timespec="seconds")})

with open(outfile, "w", newline="", encoding="utf-8-sig") as f:
    writer = csv.DictWriter(f, fieldnames=["GSE", "filename", "size_bytes", "size_MB", "md5", "download_verified_at"])
    writer.writeheader()
    writer.writerows(rows)

print("\nDownload manifest:")
for row in rows:
    print(f"  {row['filename']} | {row['size_MB']} MB | MD5: {row['md5']}")
print("\nSaved:", outfile)
PY

# 6. Inspect RAW.tar
echo ""
echo "===== 6. Inspecting GSE221561_RAW.tar ====="
python3 - <<'PY'
import tarfile, os, csv
from collections import Counter

tar_path = "01_raw/GEO/GSE221561/GSE221561_RAW.tar"
outfile = "00_metadata/03_GSE221561_RAW_archive_inventory.csv"

if not os.path.exists(tar_path):
    raise SystemExit("ERROR: GSE221561_RAW.tar was not found.")

rows = []
with tarfile.open(tar_path, "r") as tar:
    members = [m for m in tar.getmembers() if m.isfile()]
    print("Files inside RAW.tar:", len(members))
    print("\nFirst 60 files:")
    print("-" * 70)
    for i, member in enumerate(members[:60], 1):
        print(f"{i:03d} {member.name} ({member.size / 1024**2:.2f} MB)")
    
    print("\nFile types inside archive:")
    print("-" * 70)
    extensions = Counter()
    for m in members:
        name = m.name.lower()
        if name.endswith(".mtx.gz"): ext = ".mtx.gz"
        elif name.endswith(".tsv.gz"): ext = ".tsv.gz"
        elif name.endswith(".csv.gz"): ext = ".csv.gz"
        elif name.endswith(".txt.gz"): ext = ".txt.gz"
        elif name.endswith(".mtx"): ext = ".mtx"
        elif name.endswith(".csv"): ext = ".csv"
        elif name.endswith(".tsv"): ext = ".tsv"
        else: ext = os.path.splitext(name)[1] or "other"
        extensions[ext] += 1
        rows.append({"archive_path": m.name, "size_bytes": m.size, "size_MB": round(m.size / 1024**2, 3), "file_type": ext})
    
    for ext, n in sorted(extensions.items(), key=lambda x: (-x[1], x[0])):
        print(f"{ext:15s} {n}")

with open(outfile, "w", newline="", encoding="utf-8-sig") as f:
    writer = csv.DictWriter(f, fieldnames=["archive_path", "size_bytes", "size_MB", "file_type"])
    writer.writeheader()
    writer.writerows(rows)

print("\nArchive inventory saved:", outfile)
PY

# 7. Inspect cell metadata
echo ""
echo "===== 7. Inspecting cell-type metadata ====="
python3 - <<'PY'
import gzip, csv, os

path = "01_raw/GEO/GSE221561/GSE221561_metadata_celltype.csv.gz"
if not os.path.exists(path):
    raise SystemExit("ERROR: metadata_celltype.csv.gz was not found.")

with gzip.open(path, "rt", encoding="utf-8-sig", errors="replace") as f:
    reader = csv.reader(f)
    header = next(reader)
    print("Number of columns:", len(header))
    print("\nColumn names:")
    print("-" * 70)
    for i, col in enumerate(header, 1):
        print(f"{i:02d}. {col}")
    print("\nFirst 3 data rows:")
    print("-" * 70)
    for i in range(3):
        try:
            row = next(reader)
            print(f"ROW {i + 1}: {row[:min(len(row), 12)]}")
        except StopIteration:
            break
PY

# 8. Final
echo ""
echo "============================================================"
echo "STEP 7 FINISHED"
echo "============================================================"
echo ""
echo "Raw data kept untouched in: 01_raw/GEO/GSE221561/"
echo "Manifest: 00_metadata/03_GSE221561_download_manifest.csv"
echo "Archive inventory: 00_metadata/03_GSE221561_RAW_archive_inventory.csv"
echo ""
echo "IMPORTANT: RAW.tar has NOT been extracted."
